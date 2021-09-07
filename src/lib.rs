#![feature(test)]
extern crate test;

pub type Bvh2d<T> = Bvh<T, 2>;
pub type Bvh3d<T> = Bvh<T, 3>;

/// A D dimensional vector.
type Vector<const D: usize> = [f32; D];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BvhIndex(usize);

/// A ray in D dimensional space.
/// Starts at `origin` and goes in the direction `direction`.
#[derive(Clone, Copy, Debug)]
pub struct Ray<const D: usize> {
    origin: Vector<D>,
    direction: Vector<D>,
}

impl<const D: usize> Ray<D> {
    fn new(origin: Vector<D>, direction: Vector<D>) -> Self {
        Self { origin, direction }
    }
}

impl<const D: usize> From<(Vector<D>, Vector<D>)> for Ray<D> {
    fn from(tuple: (Vector<D>, Vector<D>)) -> Self {
        Self::new(tuple.0, tuple.1)
    }
}

pub trait Bounded<const D: usize>: Clone {
    type Bound: Bounded<D>;

    fn bounds(&self) -> Bounds<D>;
    fn intersect(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)>;
}

/// A bounding box in D dimensional space.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Bounds<const D: usize> {
    min: Vector<D>,
    max: Vector<D>,
}

impl<const D: usize> Bounds<D> {
    pub fn new(min: Vector<D>, max: Vector<D>) -> Bounds<D> {
        Bounds { min, max }
    }

    fn axis_length(&self, axis: usize) -> f32 {
        self.max[axis] - self.min[axis]
    }

    fn centroid(&self) -> Vector<D> {
        let mut centroid = [0.0; D];
        (0..D).for_each(|i| centroid[i] = (self.min[i] + self.max[i]) / 2.0);
        centroid
    }

    fn surface_area(&self) -> f32 {
        let mut total_area = 0.0;
        for i in 0..D {
            let length = self.axis_length(i);
            for j in (i + 1)..D {
                let length2 = self.axis_length(j);
                total_area += length * length2;
            }
        }
        (D - 1) as f32 * total_area
    }

    fn union(&self, other: &Bounds<D>) -> Bounds<D> {
        let mut min = [0.0; D];
        let mut max = [0.0; D];

        (0..D).for_each(|i| {
            min[i] = self.min[i].min(other.min[i]);
            max[i] = self.max[i].max(other.max[i])
        });

        Bounds { min: min, max: max }
    }

    fn overlap(&self, other: &Bounds<D>) -> bool {
        (0..D).all(|i| self.min[i] <= other.max[i] && self.max[i] >= other.min[i])
    }
}

impl<const D: usize> Bounded<D> for Bounds<D> {
    type Bound = Bounds<D>;

    fn bounds(&self) -> Bounds<D> {
        *self
    }

    fn intersect(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)> {
        let mut loc_t_min = t_min;
        let mut loc_t_max = t_max;

        for a in 0..D {
            let inv_d = 1.0 / ray.direction[a];
            let mut t0 = (self.min[a] - ray.origin[a]) * inv_d;
            let mut t1 = (self.max[a] - ray.origin[a]) * inv_d;

            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }
            loc_t_min = if t0 > loc_t_min { t0 } else { loc_t_min };
            loc_t_max = if t1 < loc_t_max { t1 } else { loc_t_max };
            if loc_t_max <= loc_t_min {
                return None;
            }
        }
        match (loc_t_min < 0.0, loc_t_max < 0.0) {
            (true, true) => None,
            (true, false) => Some((loc_t_max, self)),
            (false, true) => Some((loc_t_min, self)),
            (false, false) => Some((loc_t_min, self)),
        }
    }
}

/// A node in a D dimensional BVH.
#[derive(Clone, Debug)]
enum BvhNode<const D: usize> {
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

impl<T: Bounded<D>, const D: usize> Bounded<D> for Bvh<T, D> {
    type Bound = T;

    fn bounds(&self) -> Bounds<D> {
        match self.root {
            BvhNode::Node { bounds, .. } => bounds,
            BvhNode::Leaf { data } => self.objects[data.0].bounds(),
        }
    }

    fn intersect(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)> {
        let mut stack = Vec::with_capacity(32);
        stack.push(&self.root);
        let mut result = None;

        while let Some(node) = stack.pop() {
            let t_max = result.map_or(t_max, |(t, _)| t);
            match node {
                BvhNode::Node { bounds, children } => {
                    if bounds.intersect(ray, t_min, t_max).is_some() {
                        stack.push(&children[1]);
                        stack.push(&children[0]);
                    }
                }
                BvhNode::Leaf { data, .. } => {
                    let obj = &self.objects[data.0];
                    match obj.intersect(ray, t_min, t_max) {
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

impl<'a, T: Bounded<D>, const D: usize> Bvh<T, D> {
    pub fn build(objects: Vec<T>) -> Self {
        let centroids: Vec<_> = objects.iter().map(|obj| obj.bounds().centroid()).collect();
        let mut indexes: Vec<BvhIndex> = (0..objects.len()).map(|v| BvhIndex(v)).collect();

        let root = Self::_build(&objects, &centroids, &mut indexes);

        Self { objects, root }
    }

    fn _build(objects: &[T], centroids: &[Vector<D>], indexes: &mut [BvhIndex]) -> BvhNode<D> {
        let bounds = indexes
            .iter()
            .map(|i| objects[i.0].bounds())
            .reduce(|acc, b| acc.union(&b))
            .expect("No objects to build bounds.");

        match indexes.len() {
            0 => panic!("No objects given."),
            1 => BvhNode::Leaf { data: indexes[0] },
            2 => BvhNode::Node {
                bounds,
                children: [
                    Box::new(BvhNode::Leaf { data: indexes[0] }),
                    Box::new(BvhNode::Leaf { data: indexes[1] }),
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
                        let centroid1 = centroids[a.0][axis];
                        let centroid2 = centroids[b.0][axis];
                        centroid1.partial_cmp(&centroid2).unwrap()
                    });

                    for bucket in 0..NUM_BUCKETS {
                        cuts[axis][bucket] = bounds.min[axis]
                            + bounds.axis_length(axis)
                                * ((bucket + 1) as f32 / (NUM_BUCKETS + 1) as f32);

                        split_idx[axis][bucket] =
                            indexes.partition_point(|o| centroids[o.0][axis] < cuts[axis][bucket]);

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
                        let centroid1 = centroids[a.0][min_axis];
                        let centroid2 = centroids[b.0][min_axis];
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

    pub fn query_ray(&'a self, ray: &'a Ray<D>) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        let predicate =
            move |bounds: &Bounds<D>| bounds.intersect(&ray, 0.0, f32::INFINITY).is_some();
        let mut stack = Vec::with_capacity(32);
        stack.push(&self.root);

        BvhLeafIterator {
            bvh: self,
            predicate,
            stack: stack,
        }
    }

    pub fn iter_objects(&'a self) -> impl Iterator<Item = (BvhIndex, &T)> + '_ {
        let mut stack = Vec::with_capacity(32);
        stack.push(&self.root);

        BvhLeafIterator {
            bvh: self,
            predicate: |_| true,
            stack,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BvhIteratorNode<'a, T: Bounded<D>, const D: usize> {
    Node(&'a Bounds<D>),
    Leaf { id: BvhIndex, obj: &'a T },
}

pub struct BvhLeafIterator<'a, T, F, const D: usize>
where
    T: Bounded<D>,
    F: Fn(&Bounds<D>) -> bool,
{
    bvh: &'a Bvh<T, D>,
    predicate: F,
    stack: Vec<&'a BvhNode<D>>,
}

impl<'a, T, F, const D: usize> Iterator for BvhLeafIterator<'a, T, F, D>
where
    T: Bounded<D>,
    F: Fn(&Bounds<D>) -> bool,
{
    type Item = (BvhIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let node = match self.stack.pop() {
            Some(node) => node,
            None => return None,
        };

        match node {
            BvhNode::Leaf { data } => Some((*data, &self.bvh.objects[data.0])),
            BvhNode::Node { bounds, children } => {
                if (self.predicate)(bounds) {
                    self.stack.push(&children[1]);
                    self.stack.push(&children[0]);
                }
                self.next()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn bounds_axis_length_2d() {
        let b = Bounds::new([0.0, 0.0], [1.2, 1.3]);
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
    }

    #[test]
    fn bounds_axis_length_3d() {
        let b = Bounds::new([0.0, 0.0, 0.0], [1.2, 1.3, 1.4]);
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
        assert_eq!(b.axis_length(2), 1.4);
    }

    #[test]
    fn bounds_surface_area_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        assert_eq!(a.surface_area(), 1.0);
    }

    #[test]
    fn bounds_surface_area_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(a.surface_area(), 6.0);
    }

    #[test]
    fn bounds_union_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let b = Bounds::new([1.0, 1.0], [2.0, 2.0]);
        assert_eq!(a.union(&b), Bounds::new([0.0, 0.0], [2.0, 2.0]));
    }

    #[test]
    fn bounds_union_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        assert_eq!(a.union(&b), Bounds::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]));
    }

    #[test]
    fn bounds_overlap_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let b = Bounds::new([2.0, 2.0], [3.0, 3.0]);
        assert!(!a.overlap(&b));
        let c = Bounds::new([0.5, 0.5], [1.5, 1.5]);
        assert!(a.overlap(&c));
    }

    #[test]
    fn bounds_overlap_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(!a.overlap(&b));
        let c = Bounds::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        assert!(a.overlap(&c));
    }

    #[test]
    fn bounds_intersect_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let r1 = ([0.5, 2.0], [0.0, -1.0]).into(); // intersects
        let r2 = ([0.5, 0.5], [0.0, 1.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0], [1.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0], [0.0, 1.0]).into(); // miss
        assert_eq!(a.intersect(&r1, f32::MIN, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.intersect(&r2, f32::MIN, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.intersect(&r3, f32::MIN, f32::MAX), None);
        assert_eq!(a.intersect(&r4, f32::MIN, f32::MAX), None);

        let b = Bounds::new([3.0, -2.0], [4.0, -1.0]);
        let r = ([4.5, -2.5], [1.0, 1.0]).into();
        assert_eq!(b.intersect(&r, f32::MIN, f32::MAX), None);

        let b = Bounds::new([4.5, 0.0], [5.5, 1.0]);
        let r = ([4.5, -2.5], [1.0, 1.0]).into();
        assert_eq!(b.intersect(&r, f32::MIN, f32::MAX), None);
    }

    #[test]
    fn bounds_intersect_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let r1 = ([0.5, 2.0, 0.5], [0.0, -1.0, 0.0]).into(); // intersects
        let r2 = ([0.5, 0.5, 0.5], [0.0, 1.0, 0.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0, 0.5], [1.0, 0.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0, 0.5], [0.0, 1.0, 0.0]).into(); // miss
        assert_eq!(a.intersect(&r1, f32::MIN, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.intersect(&r2, f32::MIN, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.intersect(&r3, f32::MIN, f32::MAX), None);
        assert_eq!(a.intersect(&r4, f32::MIN, f32::MAX), None);
    }

    #[test]
    fn bvh_build_2d() {
        let objects = vec![
            Bounds::new([3.0, -2.0], [4.0, -1.0]),
            Bounds::new([4.5, 0.0], [5.5, 1.0]),
            Bounds::new([0.0, 1.0], [1.0, 2.0]),
            Bounds::new([1.5, -0.5], [2.5, 0.5]),
        ];
        let bvh = Bvh2d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, -2.0], [5.5, 2.0]));

        let objects = vec![
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
            Bounds::new([0.0, 0.0], [1.0, 1.0]),
        ];
        let bvh = Bvh2d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0], [1.0, 1.0]));

        let objects = vec![
            Bounds::new([3.0, -2.0], [4.0, -1.0]),
            Bounds::new([4.5, 0.0], [5.5, 1.0]),
            Bounds::new([0.0, 1.0], [1.0, 2.0]),
        ];
        let bvh = Bvh2d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, -2.0], [5.5, 2.0]));
    }

    #[test]
    fn bvh_build_3d() {
        let objects = vec![
            Bounds::new([3.0, 3.0, 3.0], [4.0, 4.0, 4.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]),
            Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
        ];
        let bvh = Bvh3d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]));

        let objects = vec![
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        ];
        let bvh = Bvh3d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]));

        let objects = vec![
            Bounds::new([3.0, 3.0, 3.0], [4.0, 4.0, 4.0]),
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]),
        ];
        let bvh = Bvh3d::build(objects);
        assert_eq!(bvh.bounds(), Bounds::new([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]));
    }

    #[test]
    fn bvh_intersect_2d() {
        let objects = vec![
            Bounds::new([3.0, -2.0], [4.0, -1.0]), // 0
            Bounds::new([4.5, 0.0], [5.5, 1.0]),   // 1
            Bounds::new([0.0, 1.0], [1.0, 2.0]),   // 2
            Bounds::new([1.5, -0.5], [2.5, 0.5]),  // 3
        ];
        let bvh = Bvh2d::build(objects);

        let r1 = ([3.5, -2.5], [0.0, 1.0]).into(); // hits 0
        let intersection = bvh.intersect(&r1, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 0.5);

        let r2 = ([6.5, 0.5], [-1.0, 0.0]).into(); // hits 1
        let intersection = bvh.intersect(&r2, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 1.0);

        let r3 = ([4.5, -2.5], [1.0, 1.0]).into(); // Hits bounding box but not a leaf
        let intersection = bvh.intersect(&r3, f32::MIN, f32::MAX);
        dbg!(&intersection);
        assert!(intersection.is_none());
    }

    #[test]
    fn bvh_intersect_3d() {
        let objects = vec![
            Bounds::new([3.0, 3.0, 3.0], [4.0, 4.0, 4.0]), // 0
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), // 1
            Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]), // 2
            Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]), // 3
        ];
        let bvh = Bvh3d::build(objects);

        let r1 = ([3.5, 5.0, 3.5], [0.0, -1.0, 0.0]).into(); // hits 0
        let intersection = bvh.intersect(&r1, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 1.0);

        let r2 = ([0.5, 0.5, -0.5], [0.0, 0.0, 1.0]).into(); // hits 1
        let intersection = bvh.intersect(&r2, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 0.5);

        let r3 = ([0.5, -0.5, 3.5], [0.0, 1.0, 0.0]).into(); // Hits bounding box but not a leaf
        let intersection = bvh.intersect(&r3, f32::MIN, f32::MAX);
        assert!(intersection.is_none());
    }

    #[test]
    fn bvh_query_ray_2d() {
        let objects = vec![
            Bounds::new([3.0, -2.0], [4.0, -1.0]), // 0
            Bounds::new([4.5, 0.0], [5.5, 1.0]),   // 1
            Bounds::new([0.0, 1.0], [1.0, 2.0]),   // 2
            Bounds::new([1.5, -0.5], [2.5, 0.5]),  // 3
        ];
        let bvh = Bvh2d::build(objects);

        let r = Ray::new([0.0, 0.0], [1.0, 0.1]);
        let mut ray_hits = bvh.query_ray(&r);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(2), &bvh.objects[2]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(3), &bvh.objects[3]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert!(ray_hits.next().is_none());
    }

    #[test]
    fn bvh_query_ray_3d() {
        let objects = vec![
            Bounds::new([3.0, 3.0, 3.0], [4.0, 4.0, 4.0]), // 0
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), // 1
            Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]), // 2
            Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]), // 3
        ];
        let bvh = Bvh3d::build(objects);

        let r = Ray::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let mut ray_hits = bvh.query_ray(&r);
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(1), &bvh.objects[1]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(3), &bvh.objects[3]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(2), &bvh.objects[2]));
        assert_eq!(ray_hits.next().unwrap(), (BvhIndex(0), &bvh.objects[0]));
        assert!(ray_hits.next().is_none());
    }

    #[bench]
    fn bench_bvh_build_2d(b: &mut Bencher) {
        let objects: Vec<Bounds<2>> = (0..100)
            .map(|i| Bounds::new([i as f32; 2], [(i + 1) as f32; 2]))
            .collect();
        b.iter(|| Bvh2d::build(objects.clone()));
    }

    #[bench]
    fn bench_bvh_build_3d(b: &mut Bencher) {
        let objects: Vec<Bounds<3>> = (0..100)
            .map(|i| Bounds::new([i as f32; 3], [(i + 1) as f32; 3]))
            .collect();
        b.iter(|| Bvh3d::build(objects.clone()));
    }

    #[bench]
    fn bench_bvh_intersect_2d(b: &mut Bencher) {
        let objects: Vec<Bounds<2>> = (0..100)
            .map(|i| Bounds::new([i as f32; 2], [(i + 1) as f32; 2]))
            .collect();
        let bvh = Bvh2d::build(objects);
        let r = ([-0.5, 49.5], [1.0, 0.0]).into();
        b.iter(|| {
            assert_eq!(
                bvh.intersect(&r, f32::MIN, f32::MAX).unwrap().1,
                &Bounds::new([49.0; 2], [50.0; 2])
            )
        });
    }

    #[bench]
    fn bench_bvh_intersect_3d(b: &mut Bencher) {
        let objects: Vec<Bounds<3>> = (0..100)
            .map(|i| Bounds::new([i as f32; 3], [(i + 1) as f32; 3]))
            .collect();
        let bvh = Bvh3d::build(objects);
        let r = ([-0.5, 49.5, 49.5], [1.0, 0.0, 0.0]).into();
        b.iter(|| {
            assert_eq!(
                bvh.intersect(&r, f32::MIN, f32::MAX).unwrap().1,
                &Bounds::new([49.0; 3], [50.0; 3])
            )
        });
    }

    #[bench]
    fn bench_bvh_iter_2d(b: &mut Bencher) {
        let objects: Vec<Bounds<2>> = (0..100)
            .map(|i| Bounds::new([i as f32; 2], [(i + 1) as f32; 2]))
            .collect();
        let bvh = Bvh2d::build(objects);
        b.iter(|| bvh.iter_objects().map(|_| ()).count());
    }

    #[bench]
    fn bench_bvh_iter_3d(b: &mut Bencher) {
        let objects: Vec<Bounds<3>> = (0..100)
            .map(|i| Bounds::new([i as f32; 3], [(i + 1) as f32; 3]))
            .collect();
        let bvh = Bvh3d::build(objects);
        b.iter(|| bvh.iter_objects().map(|_| ()).count());
    }

    #[bench]
    fn bench_bvh_query_ray_2d(b: &mut Bencher) {
        let objects: Vec<Bounds<2>> = (0..100)
            .map(|i| Bounds::new([i as f32; 2], [(i + 1) as f32; 2]))
            .collect();
        let bvh = Bvh2d::build(objects);
        let r = ([-0.5, 49.5], [1.0, 0.0]).into();
        b.iter(|| assert!(bvh.query_ray(&r).any(|(idx, _)| idx == BvhIndex(49))));
    }

    #[bench]
    fn bench_bvh_query_ray_3d(b: &mut Bencher) {
        let objects: Vec<Bounds<3>> = (0..100)
            .map(|i| Bounds::new([i as f32; 3], [(i + 1) as f32; 3]))
            .collect();
        let bvh = Bvh3d::build(objects);
        let r = ([-0.5, 49.5, 49.5], [1.0, 0.0, 0.0]).into();
        b.iter(|| assert!(bvh.query_ray(&r).any(|(idx, _)| idx == BvhIndex(49))));
    }
}
