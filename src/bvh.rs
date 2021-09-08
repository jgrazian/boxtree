use std::ops::Index;

use crate::bounds::Bounds;
use crate::iter::BvhLeafIterator;
use crate::traits::*;
use crate::{Ray, Vector};

pub type Bvh2d<T> = Bvh<T, 2>;
pub type Bvh3d<T> = Bvh<T, 3>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BvhIndex(usize);

/// A node in a D dimensional BVH.
#[derive(Clone, Debug)]
pub(crate) enum BvhNode<const D: usize> {
    Node {
        bounds: Bounds<D>,
        children: [Box<BvhNode<D>>; 2],
    },
    Leaf {
        data: BvhIndex,
    },
}

/// A D dimensional bounding volume hierarchy.
#[derive(Clone, Debug)]
pub struct Bvh<T: Bounded<D>, const D: usize> {
    objects: Vec<T>,
    root: BvhNode<D>,
}

impl<'a, T: Bounded<D>, const D: usize> Bvh<T, D> {
    const STACK_SIZE: usize = 32;

    pub fn build(objects: Vec<T>) -> Self {
        let centroids: Vec<_> = objects.iter().map(|obj| obj.bounds().centroid()).collect();
        let mut indexes: Vec<usize> = (0..objects.len()).collect();

        let root = Self::_build(&objects, &centroids, &mut indexes);

        Self { objects, root }
    }

    fn _build(objects: &[T], centroids: &[Vector<D>], indexes: &mut [usize]) -> BvhNode<D> {
        let bounds = indexes
            .iter()
            .map(|i| objects[*i].bounds())
            .reduce(|acc, b| acc.union(&b))
            .expect("No objects to build bounds.");

        match indexes.len() {
            0 => panic!("No objects given."),
            1 => BvhNode::Leaf {
                data: BvhIndex(indexes[0]),
            },
            2 => BvhNode::Node {
                bounds,
                children: [
                    Box::new(BvhNode::Leaf {
                        data: BvhIndex(indexes[0]),
                    }),
                    Box::new(BvhNode::Leaf {
                        data: BvhIndex(indexes[1]),
                    }),
                ],
            },
            _ => {
                const NUM_BUCKETS: usize = 16;
                let mut cuts = [[0.0; NUM_BUCKETS]; D];
                let mut split_idx = [[0 as usize; NUM_BUCKETS]; D];
                let mut scores = [[0.0; NUM_BUCKETS]; D];
                let bounds_sa = bounds.surface_area();

                for axis in 0..D {
                    // Sort objects by axis
                    indexes.sort_unstable_by(|a, b| {
                        let centroid1 = centroids[*a][axis];
                        let centroid2 = centroids[*b][axis];
                        centroid1.partial_cmp(&centroid2).unwrap()
                    });

                    for bucket in 0..NUM_BUCKETS {
                        cuts[axis][bucket] = bounds.min[axis]
                            + bounds.axis_length(axis)
                                * ((bucket + 1) as f32 / (NUM_BUCKETS + 1) as f32);

                        split_idx[axis][bucket] =
                            indexes.partition_point(|o| centroids[*o][axis] < cuts[axis][bucket]);

                        let (left, right) = indexes.split_at(split_idx[axis][bucket]);

                        // Bad cut location one of the sides has no shapes.
                        if left.len() == 0 || right.len() == 0 {
                            scores[axis][bucket] = f32::INFINITY;
                            continue;
                        }

                        let mut left_bounds = bounds.clone();
                        left_bounds.max[axis] = cuts[axis][bucket];
                        let mut right_bounds = bounds.clone();
                        right_bounds.min[axis] = cuts[axis][bucket];

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
                for axis in 0..D {
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
                    indexes.len() / 2 // Centroids are all very close just split in half
                } else {
                    indexes.sort_unstable_by(|a, b| {
                        let centroid1 = centroids[*a][min_axis];
                        let centroid2 = centroids[*b][min_axis];
                        centroid1.partial_cmp(&centroid2).unwrap()
                    });

                    split_idx[min_axis][min_bucket]
                };

                let (mut left, mut right) = indexes.split_at_mut(split_index);

                let left_node = Self::_build(objects, centroids, &mut left);
                let right_node = Self::_build(objects, centroids, &mut right);

                BvhNode::Node {
                    bounds,
                    children: [Box::new(left_node), Box::new(right_node)],
                }
            }
        }
    }

    pub fn query_ray(
        &'a self,
        ray: &'a Ray<D>,
        t_min: f32,
        t_max: f32,
    ) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        let predicate = move |bounds: &Bounds<D>| bounds.ray_hit(&ray, t_min, t_max).is_some();

        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.root);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack,
        }
    }

    pub fn query_bounds(
        &'a self,
        bounds: &'a Bounds<D>,
    ) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        let predicate = move |bounds2: &Bounds<D>| bounds.bounds_hit(&bounds2);

        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.root);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack,
        }
    }

    pub fn query_point(
        &'a self,
        point: &'a Vector<D>,
    ) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        let predicate = move |bounds: &Bounds<D>| bounds.point_hit(&point);

        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.root);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack,
        }
    }

    pub fn iter_objects(&'a self) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.root);

        BvhLeafIterator {
            bvh: self,
            predicate: |_| true,
            stack,
        }
    }
}

impl<'a, T: RayHittable<D>, const D: usize> Bvh<T, D> {
    pub fn query_ray_exact(
        &'a self,
        ray: &'a Ray<D>,
        t_min: f32,
        t_max: f32,
    ) -> impl Iterator<Item = (f32, BvhIndex, &T)> + '_ {
        self.query_ray(ray, t_min, t_max)
            .filter_map(move |(idx, obj)| match obj.ray_hit(ray, t_min, t_max) {
                None => None,
                Some((t, _)) => Some((t, idx, obj)),
            })
    }
}

impl<'a, T: BoundsHittable<D>, const D: usize> Bvh<T, D> {
    pub fn query_bounds_exact(
        &'a self,
        bounds: &'a Bounds<D>,
    ) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        self.query_bounds(bounds)
            .filter(move |(_, obj)| obj.bounds_hit(bounds))
    }
}

impl<'a, T: PointHittable<D>, const D: usize> Bvh<T, D> {
    pub fn query_point_exact(
        &'a self,
        point: &'a Vector<D>,
    ) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        self.query_point(point)
            .filter(move |(_, obj)| obj.point_hit(point))
    }
}

impl<T: Bounded<D>, const D: usize> Index<BvhIndex> for Bvh<T, D> {
    type Output = T;

    fn index(&self, index: BvhIndex) -> &Self::Output {
        &self.objects[index.0]
    }
}

impl<T: Bounded<D>, const D: usize> Index<&BvhIndex> for Bvh<T, D> {
    type Output = T;

    fn index(&self, index: &BvhIndex) -> &Self::Output {
        &self.objects[index.0]
    }
}

impl<T: Bounded<D>, const D: usize> Bounded<D> for Bvh<T, D> {
    type Bound = T;

    fn bounds(&self) -> Bounds<D> {
        match self.root {
            BvhNode::Node { bounds, .. } => bounds,
            BvhNode::Leaf { data } => self[data].bounds(),
        }
    }
}

impl<T: RayHittable<D>, const D: usize> RayHittable<D> for Bvh<T, D> {
    fn ray_hit(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &T)> {
        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.root);
        let mut result = None;

        while let Some(node) = stack.pop() {
            let t_max = result.map_or(t_max, |(t, _)| t);
            match node {
                BvhNode::Node { bounds, children } => {
                    if bounds.ray_hit(ray, t_min, t_max).is_some() {
                        stack.push(&children[1]);
                        stack.push(&children[0]);
                    }
                }
                BvhNode::Leaf { data, .. } => {
                    let obj = &self[data];
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
    /// ┼──────────
    fn setup_2d() -> Bvh2d<Bounds<2>> {
        let objects = vec![
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([1.0, 1.0], [2.0, 2.0]),
            Bounds::new([3.0, 0.0], [4.0, 1.0]),
            Bounds::new([4.0, 1.0], [5.0, 2.0]),
        ];
        Bvh2d::build(objects)
    }

    fn setup_3d() -> Bvh3d<Bounds<3>> {
        let objects = vec![
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
            Bounds::new([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]),
            Bounds::new([4.0, 1.0, 1.0], [5.0, 2.0, 2.0]),
        ];
        Bvh3d::build(objects)
    }

    #[test]
    fn build_2d() {
        let bvh = setup_2d();
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0], [5.0, 2.0]));

        let objects = vec![
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
        ];
        let bvh = Bvh2d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0], [1.0, 1.0]));
    }

    #[test]
    fn build_3d() {
        let bvh = setup_3d();
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0, 0.0], [5.0, 2.0, 2.0]));

        let objects = vec![
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        ];
        let bvh = Bvh3d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]));
    }

    #[test]
    fn query_ray_2d() {
        let bvh = setup_2d();

        let r = Ray::new([1.0, 0.0], [0.0, 1.0]);
        let mut ray_hits = bvh.query_ray(&r, 0.0, f32::INFINITY);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_ray_3d() {
        let bvh = setup_3d();

        let r = Ray::new([1.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let mut ray_hits = bvh.query_ray(&r, 0.0, f32::INFINITY);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_ray_exact_2d() {
        let bvh = setup_2d();

        let r = Ray::new([0.0, 0.5], [1.0, 0.0]);
        let mut ray_hits = bvh.query_ray_exact(&r, 0.0, f32::MAX);
        assert_eq!(
            ray_hits.next().unwrap(),
            (0.0, BvhIndex(0), &bvh.objects[0])
        );
        assert_eq!(
            ray_hits.next().unwrap(),
            (3.0, BvhIndex(2), &bvh.objects[2])
        );
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_ray_exact_3d() {
        let bvh = setup_3d();

        let r = Ray::new([0.0, 0.5, 0.5], [1.0, 0.0, 0.0]);
        let mut ray_hits = bvh.query_ray_exact(&r, 0.0, f32::MAX);
        assert_eq!(
            ray_hits.next().unwrap(),
            (0.0, BvhIndex(0), &bvh.objects[0])
        );
        assert_eq!(
            ray_hits.next().unwrap(),
            (3.0, BvhIndex(2), &bvh.objects[2])
        );
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_bounds_2d() {
        let bvh = setup_2d();

        let b = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let mut ray_hits = bvh.query_bounds(&b);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_bounds_3d() {
        let bvh = setup_3d();

        let b = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mut ray_hits = bvh.query_bounds(&b);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_bounds_exact_2d() {
        let bvh = setup_2d();

        let b = Bounds::new([0.0, 0.0], [5.0, 0.5]);
        let mut ray_hits = bvh.query_bounds_exact(&b);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(2), &bvh.objects[2]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_bounds_exact_3d() {
        let bvh = setup_3d();

        let b = Bounds::new([0.0, 0.0, 0.0], [5.0, 0.5, 0.5]);
        let mut ray_hits = bvh.query_bounds_exact(&b);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(2), &bvh.objects[2]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_point_2d() {
        let bvh = setup_2d();

        let b = [0.5, 0.5];
        let mut ray_hits = bvh.query_point(&b);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_point_3d() {
        let bvh = setup_3d();

        let b = [0.5, 0.5, 0.5];
        let mut ray_hits = bvh.query_point(&b);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_point_exact_2d() {
        let bvh = setup_2d();

        let p = [0.5, 0.5];
        let mut ray_hits = bvh.query_point_exact(&p);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn query_point_exact_3d() {
        let bvh = setup_3d();

        let p = [0.5, 0.5, 0.5];
        let mut ray_hits = bvh.query_point_exact(&p);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert!(ray_hits.next().is_none());
    }
}
