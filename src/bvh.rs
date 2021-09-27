use std::ops::Index;

use crate::bounds::{Bounds2, Bounds3, Bounds3A};
use crate::iter::BvhLeafIterator;
use crate::traits::*;

pub type Bvh2<T> = Bvh<Bounds2, T, 2>;
pub type Bvh3<T> = Bvh<Bounds3, T, 2>;
pub type Bvh3A<T> = Bvh<Bounds3A, T, 2>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BvhObjKey(usize);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BvhNodeKey(pub(crate) usize);

/// A node in a D dimensional BVH.
#[derive(Clone, Debug)]
pub(crate) enum BvhNode<B: BoundingBox, const N: usize> {
    Node {
        bounds: B,
        children: [Option<BvhNodeKey>; N],
    },
    Leaf(BvhObjKey),
}

/// A D dimensional bounding volume hierarchy.
#[derive(Clone, Debug)]
pub struct Bvh<B: BoundingBox, T: Bounded<B>, const N: usize> {
    objects: Vec<T>,
    pub(crate) nodes: Vec<BvhNode<B, N>>,
}

impl<'a, B: BoundingBox, T: Bounded<B>, const N: usize> Bvh<B, T, N> {
    const STACK_SIZE: usize = 32;

    pub fn build(mut objects: Vec<T>) -> Self {
        let centroids: Vec<_> = objects.iter().map(|obj| obj.bounds().centroid()).collect();
        let mut indexes = (0..objects.len()).collect::<Vec<_>>();

        let mut nodes = Vec::with_capacity(Self::STACK_SIZE);
        let mut stack: Vec<(Option<(BvhNodeKey, usize)>, &mut [usize])> =
            vec![(None, &mut indexes)];

        while let Some((parent, idxs)) = stack.pop() {
            match idxs.len() {
                0 => (),
                1 => match parent {
                    None => nodes.push(BvhNode::Leaf(BvhObjKey(idxs[0]))),
                    Some((parent_key, pos)) => {
                        nodes.push(BvhNode::Leaf(BvhObjKey(idxs[0])));
                        let leaf_id = BvhNodeKey(nodes.len() - 1);

                        match nodes[parent_key.0] {
                            BvhNode::Node {
                                ref mut children, ..
                            } => children[pos] = Some(leaf_id),
                            _ => panic!("Expected node found leaf."),
                        }
                    }
                },
                _ => {
                    let bounds = idxs
                        .iter()
                        .map(|i| objects[*i].bounds())
                        .reduce(|acc, b| acc.union(&b))
                        .expect("No objects to build bounds.");

                    nodes.push(BvhNode::Node {
                        bounds,
                        children: [None; N],
                    });
                    let node_id = BvhNodeKey(nodes.len() - 1);

                    match parent {
                        None => (),
                        Some((parent_key, pos)) => match nodes[parent_key.0] {
                            BvhNode::Node {
                                ref mut children, ..
                            } => children[pos] = Some(node_id),
                            _ => panic!("Expected node found leaf."),
                        },
                    }

                    let mut chunks = if N == 2 {
                        Self::_split_sah(&bounds, &objects, &centroids, idxs)
                    } else {
                        Self::_split_chunks(&bounds, &centroids, idxs)
                    };

                    let num_chunks = chunks.iter().filter(|c| c.is_some()).count();
                    for i in (0..num_chunks).rev() {
                        let index_slice = std::mem::take(&mut chunks[i]).unwrap();
                        stack.push((Some((node_id, i)), index_slice));
                    }
                }
            }
        }

        let mut sorted_obj_ids = nodes
            .iter()
            .filter_map(|n| match n {
                BvhNode::Node { .. } => None,
                BvhNode::Leaf(obj_key) => Some(obj_key),
            })
            .enumerate()
            .collect::<Vec<_>>();
        sorted_obj_ids.sort_by_key(|v| v.1);

        let mut tmp_objs = sorted_obj_ids
            .iter()
            .map(|v| v.0)
            .zip(objects)
            .collect::<Vec<_>>();
        tmp_objs.sort_by_key(|v| v.0);

        objects = tmp_objs.into_iter().map(|v| v.1).collect::<Vec<_>>();
        nodes
            .iter_mut()
            .filter_map(|n| match n {
                BvhNode::Node { .. } => None,
                BvhNode::Leaf(ref mut obj_key) => Some(obj_key),
            })
            .enumerate()
            .for_each(|(i, leaf)| leaf.0 = i);

        Self { objects, nodes }
    }

    #[inline]
    fn _split_sah(
        bounds: &B,
        objects: &[T],
        centroids: &'a [B::Vector],
        indexes: &'a mut [usize],
    ) -> [Option<&'a mut [usize]>; N] {
        const MAX_DIM: usize = 3;
        const NUM_BUCKETS: usize = 16;
        let mut cuts = [[0.0; NUM_BUCKETS]; MAX_DIM];
        let mut split_idx = [[0 as usize; NUM_BUCKETS]; MAX_DIM];
        let mut scores = [[f32::INFINITY; NUM_BUCKETS]; MAX_DIM];
        let bounds_sa = bounds.surface_area();

        for axis in 0..B::DIM {
            // Sort objects by axis
            indexes.sort_unstable_by(|a, b| {
                let centroid1 = centroids[*a][axis];
                let centroid2 = centroids[*b][axis];
                centroid1.partial_cmp(&centroid2).unwrap()
            });

            for bucket in 0..NUM_BUCKETS {
                cuts[axis][bucket] = bounds.min()[axis]
                    + bounds.axis_length(axis) * ((bucket + 1) as f32 / (NUM_BUCKETS + 1) as f32);

                split_idx[axis][bucket] =
                    indexes.partition_point(|o| centroids[*o][axis] < cuts[axis][bucket]);

                let (left, right) = indexes.split_at(split_idx[axis][bucket]);

                // Bad cut location one of the sides has no shapes.
                if left.len() == 0 || right.len() == 0 {
                    continue;
                }

                let mut left_bounds = bounds.clone();
                left_bounds.max_mut()[axis] = objects[left[left.len() - 1]].bounds().max()[axis];
                let mut right_bounds = bounds.clone();
                right_bounds.min_mut()[axis] = objects[right[0]].bounds().min()[axis];

                let left_sa = left_bounds.surface_area();
                let right_sa = right_bounds.surface_area();

                let left_ratio = left_sa / bounds_sa;
                let right_ratio = right_sa / bounds_sa;

                scores[axis][bucket] =
                    1.0 + left_ratio * left.len() as f32 + right_ratio * right.len() as f32;
            }
        }

        // Select minimum score
        let mut min_score = f32::INFINITY;
        let mut min_axis = 0;
        let mut min_bucket = 0;
        for axis in 0..B::DIM {
            for bucket in 0..NUM_BUCKETS {
                if scores[axis][bucket] < min_score {
                    min_score = scores[axis][bucket];
                    min_axis = axis;
                    min_bucket = bucket;
                }
            }
        }

        // Decide where to split
        let split_index = if min_score == f32::INFINITY {
            // Centroids are all very close just split in half
            indexes.len() / 2
        } else {
            indexes.sort_unstable_by(|a, b| {
                let centroid1 = centroids[*a][min_axis];
                let centroid2 = centroids[*b][min_axis];
                centroid1.partial_cmp(&centroid2).unwrap()
            });

            split_idx[min_axis][min_bucket]
        };

        let (left, right) = indexes.split_at_mut(split_index);

        #[allow(dead_code)]
        const TEMP: Option<&mut [usize]> = None;
        let mut out = [TEMP; N];
        out[0] = Some(left);
        out[1] = Some(right);
        out
    }

    #[inline]
    fn _split_chunks(
        bounds: &B,
        centroids: &'a [B::Vector],
        indexes: &'a mut [usize],
    ) -> [Option<&'a mut [usize]>; N] {
        let mut min = (0, f32::MIN);
        for i in 0..B::DIM {
            let side_len = bounds.axis_length(i);
            if side_len > min.1 {
                min = (i, side_len);
            }
        }
        let longest_axis_idx = min.0;

        indexes.sort_unstable_by(|a, b| {
            let centroid1 = centroids[*a][longest_axis_idx];
            let centroid2 = centroids[*b][longest_axis_idx];
            centroid1.partial_cmp(&centroid2).unwrap()
        });

        let chunk_size = match indexes.len() {
            n if n <= N => 1,
            n if n <= (N * N) => N,
            n => (n + (N - 1)) / N,
        };

        #[allow(dead_code)]
        const TEMP: Option<&mut [usize]> = None;
        let mut out = [TEMP; N];
        indexes
            .chunks_mut(chunk_size)
            .zip(out.iter_mut())
            .for_each(|(chunk, out_chunk)| *out_chunk = Some(chunk));
        out
    }

    pub fn query_ray(
        &'a self,
        ray: &'a B::Ray,
        t_min: f32,
        t_max: f32,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let predicate = move |bounds: &B| bounds.ray_hit(&ray, t_min, t_max).is_some();

        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.nodes[0]);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack,
        }
    }

    pub fn query_bounds(&'a self, bounds: &'a B) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let predicate = move |bounds2: &B| bounds.bounds_hit(&bounds2);

        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.nodes[0]);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack,
        }
    }

    pub fn query_point(
        &'a self,
        point: &'a B::Vector,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let predicate = move |bounds: &B| bounds.point_hit(&point);

        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.nodes[0]);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack,
        }
    }

    pub fn iter_objects(&'a self) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.nodes[0]);

        BvhLeafIterator {
            bvh: self,
            predicate: |_| true,
            stack,
        }
    }
}

impl<'a, B: BoundingBox, T: RayHittable<B>, const N: usize> Bvh<B, T, N> {
    pub fn query_ray_exact(
        &'a self,
        ray: &'a B::Ray,
        t_min: f32,
        t_max: f32,
    ) -> impl Iterator<Item = (f32, BvhObjKey, &T)> + '_ {
        self.query_ray(ray, t_min, t_max)
            .filter_map(move |(idx, obj)| match obj.ray_hit(ray, t_min, t_max) {
                None => None,
                Some((t, _)) => Some((t, idx, obj)),
            })
    }
}

impl<'a, B: BoundingBox, T: BoundsHittable<B>, const N: usize> Bvh<B, T, N> {
    pub fn query_bounds_exact(
        &'a self,
        bounds: &'a B,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        self.query_bounds(bounds)
            .filter(move |(_, obj)| obj.bounds_hit(bounds))
    }
}

impl<'a, B: BoundingBox, T: PointHittable<B>, const N: usize> Bvh<B, T, N> {
    pub fn query_point_exact(
        &'a self,
        point: &'a B::Vector,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        self.query_point(point)
            .filter(move |(_, obj)| obj.point_hit(point))
    }
}

impl<B: BoundingBox, T: Bounded<B>, const N: usize> Index<BvhObjKey> for Bvh<B, T, N> {
    type Output = T;

    fn index(&self, index: BvhObjKey) -> &Self::Output {
        &self.objects[index.0]
    }
}

impl<B: BoundingBox, T: Bounded<B>, const N: usize> Index<&BvhObjKey> for Bvh<B, T, N> {
    type Output = T;

    fn index(&self, index: &BvhObjKey) -> &Self::Output {
        &self.objects[index.0]
    }
}

impl<B: BoundingBox, T: Bounded<B>, const N: usize> Bounded<B> for Bvh<B, T, N> {
    type Item = T;

    fn bounds(&self) -> B {
        match self.nodes[0] {
            BvhNode::Node { bounds, .. } => bounds,
            BvhNode::Leaf(obj_key) => self[obj_key].bounds(),
        }
    }
}

impl<B: BoundingBox, T: RayHittable<B>, const N: usize> RayHittable<B> for Bvh<B, T, N> {
    fn ray_hit(&self, ray: &B::Ray, t_min: f32, t_max: f32) -> Option<(f32, &T)> {
        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.nodes[0]);
        let mut result = None;

        while let Some(node) = stack.pop() {
            let t_max = result.map_or(t_max, |(t, _)| t);
            match node {
                BvhNode::Node { bounds, children } => {
                    if bounds.ray_hit(ray, t_min, t_max).is_some() {
                        for child_key in children.iter().rev() {
                            match child_key {
                                Some(key) => stack.push(&self.nodes[key.0]),
                                None => (),
                            }
                        }
                    }
                }
                BvhNode::Leaf(obj_key) => {
                    let obj = &self[obj_key];
                    match obj.ray_hit(ray, t_min, t_max) {
                        None => (),
                        Some((t, _)) => {
                            result = Some((t, obj));
                        }
                    }
                }
            }
        }
        result
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// │╔══╤═╗ ╔══╤═╗
    /// │╟─┐└─╢ ╟─┐└─╢
    /// │╚═╧══╝ ╚═╧══╝
    /// ┼───────────────
    fn setup_2d() -> Bvh2<Bounds2> {
        let objects = vec![
            ([3.0, 0.0], [4.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([4.0, 1.0], [5.0, 2.0]).into(),
            ([1.0, 1.0], [2.0, 2.0]).into(),
        ];
        Bvh2::build(objects)
    }

    fn setup_3d() -> Bvh3<Bounds3> {
        let objects = vec![
            ([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([4.0, 1.0, 1.0], [5.0, 2.0, 2.0]).into(),
            ([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]).into(),
        ];
        Bvh3::build(objects)
    }

    #[test]
    fn build_2d() {
        let bvh = setup_2d();
        assert_eq!(bvh.bounds(), ([0.0, 0.0], [5.0, 2.0]).into());

        let objects: Vec<Bounds2> = vec![
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
        ];
        let bvh = Bvh2::build(objects);
        assert_eq!(bvh.bounds(), ([0.0, 0.0], [1.0, 1.0]).into());
    }

    #[test]
    fn build_3d() {
        let bvh = setup_3d();
        assert_eq!(bvh.bounds(), ([0.0, 0.0, 0.0], [5.0, 2.0, 2.0]).into());

        let objects: Vec<Bounds3> = vec![
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
        ];
        let bvh = Bvh3::build(objects);
        assert_eq!(bvh.bounds(), ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into());
    }

    #[test]
    fn query_ray_2d() {
        let bvh = setup_2d();

        let r = ([1.0, 0.0], [0.0, 1.0]).into();
        let mut ray_hits = bvh.query_ray(&r, 0.0, f32::INFINITY);
        assert_eq!(ray_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhObjKey(1), &bvh.objects[1]));
        assert_eq!(ray_hits.next(), None);
    }

    #[test]
    fn query_ray_3d() {
        let bvh = setup_3d();

        let r = ([1.0, 1.0, 0.0], [0.0, 0.0, 1.0]).into();
        let mut ray_hits = bvh.query_ray(&r, 0.0, f32::INFINITY);
        assert_eq!(ray_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhObjKey(1), &bvh.objects[1]));
        assert_eq!(ray_hits.next(), None);
    }

    #[test]
    fn query_ray_exact_2d() {
        let bvh = setup_2d();

        let r = ([0.0, 0.5], [1.0, 0.0]).into();
        let mut ray_hits = bvh.query_ray_exact(&r, 0.0, f32::MAX);
        assert_eq!(
            ray_hits.next().unwrap(),
            (0.0, BvhObjKey(0), &bvh.objects[0])
        );
        assert_eq!(
            ray_hits.next().unwrap(),
            (3.0, BvhObjKey(2), &bvh.objects[2])
        );
        assert_eq!(ray_hits.next(), None);
    }

    #[test]
    fn query_ray_exact_3d() {
        let bvh = setup_3d();

        let r = ([0.0, 0.5, 0.5], [1.0, 0.0, 0.0]).into();
        let mut ray_hits = bvh.query_ray_exact(&r, 0.0, f32::MAX);
        assert_eq!(
            ray_hits.next().unwrap(),
            (0.0, BvhObjKey(0), &bvh.objects[0])
        );
        assert_eq!(
            ray_hits.next().unwrap(),
            (3.0, BvhObjKey(2), &bvh.objects[2])
        );
        assert_eq!(ray_hits.next(), None);
    }

    #[test]
    fn query_bounds_2d() {
        let bvh = setup_2d();

        let b = ([0.0, 0.0], [1.0, 1.0]).into();
        let mut bounds_hits = bvh.query_bounds(&b);
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(1), &bvh.objects[1]));
        assert_eq!(bounds_hits.next(), None);
    }

    #[test]
    fn query_bounds_3d() {
        let bvh = setup_3d();

        let b = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let mut bounds_hits = bvh.query_bounds(&b);
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(1), &bvh.objects[1]));
        assert_eq!(bounds_hits.next(), None);
    }

    #[test]
    fn query_bounds_exact_2d() {
        let bvh = setup_2d();

        let b = ([0.0, 0.0], [5.0, 0.5]).into();
        let mut bounds_hits = bvh.query_bounds_exact(&b);
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(2), &bvh.objects[2]));
        assert_eq!(bounds_hits.next(), None);
    }

    #[test]
    fn query_bounds_exact_3d() {
        let bvh = setup_3d();

        let b = ([0.0, 0.0, 0.0], [5.0, 0.5, 0.5]).into();
        let mut bounds_hits = bvh.query_bounds_exact(&b);
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(2), &bvh.objects[2]));
        assert_eq!(bounds_hits.next(), None);
    }

    #[test]
    fn query_point_2d() {
        let bvh = setup_2d();

        let b = [0.5, 0.5].into();
        let mut point_hits = bvh.query_point(&b);
        assert_eq!(point_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(point_hits.next().unwrap(), (BvhObjKey(1), &bvh.objects[1]));
        assert_eq!(point_hits.next(), None);
    }

    #[test]
    fn query_point_3d() {
        let bvh = setup_3d();

        let b = [0.5, 0.5, 0.5].into();
        let mut point_hits = bvh.query_point(&b);
        assert_eq!(point_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(point_hits.next().unwrap(), (BvhObjKey(1), &bvh.objects[1]));
        assert_eq!(point_hits.next(), None);
    }

    #[test]
    fn query_point_exact_2d() {
        let bvh = setup_2d();

        let p = [0.5, 0.5].into();
        let mut point_hits = bvh.query_point_exact(&p);
        assert_eq!(point_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(point_hits.next(), None);
    }

    #[test]
    fn query_point_exact_3d() {
        let bvh = setup_3d();

        let p = [0.5, 0.5, 0.5].into();
        let mut point_hits = bvh.query_point_exact(&p);
        assert_eq!(point_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(point_hits.next(), None);
    }
}
