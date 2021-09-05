#![feature(test)]
extern crate test;

pub type Bvh2d<T> = Bvh<T, 2>;
pub type Bvh3d<T> = Bvh<T, 3>;

type Vector<const D: usize> = [f32; D];
type Ray<const D: usize> = (Vector<D>, Vector<D>);

pub trait Bounded<const D: usize>: Clone {
    type Bound: Bounded<D>;

    fn bounds(&self) -> Bounds<D>;
    fn intersect(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)>;
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Bounds<const D: usize> {
    min: Vector<D>,
    max: Vector<D>,
}

impl<const D: usize> Bounds<D> {
    pub fn new(min: Vector<D>, max: Vector<D>) -> Bounds<D> {
        Bounds { min, max }
    }

    fn length(&self, axis: usize) -> f32 {
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
            let length = self.length(i);
            for j in (i + 1)..D {
                let length2 = self.length(j);
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
            let inv_d = 1.0 / ray.1[a];
            let mut t0 = (self.min[a] - ray.0[a]) * inv_d;
            let mut t1 = (self.max[a] - ray.0[a]) * inv_d;

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

#[derive(Clone, Debug)]
enum BvhNode<const D: usize> {
    Node {
        bounds: Bounds<D>,
        children: [Box<BvhNode<D>>; 2],
    },
    Leaf {
        bounds: Bounds<D>,
        data: usize,
    },
}

impl<const D: usize> BvhNode<D> {
    fn intersect(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, usize)> {
        match self {
            BvhNode::Node { bounds, children } => match bounds.intersect(ray, t_min, t_max) {
                None => None,
                Some(_) => {
                    let left_hit = children[0].intersect(ray, t_min, t_max);
                    let right_hit = children[1].intersect(ray, t_min, t_max);
                    match (left_hit, right_hit) {
                        (None, None) => None,
                        (Some((t, idx)), None) => Some((t, idx)),
                        (None, Some((t, idx))) => Some((t, idx)),
                        (Some((t0, idx0)), Some((t1, idx1))) => {
                            if t0 < t1 {
                                Some((t0, idx0))
                            } else {
                                Some((t1, idx1))
                            }
                        }
                    }
                }
            },
            BvhNode::Leaf { bounds, data } => match bounds.intersect(ray, t_min, t_max) {
                None => None,
                Some((t, _)) => Some((t, *data)),
            },
        }
    }
}

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
            BvhNode::Leaf { bounds, .. } => bounds,
        }
    }

    fn intersect(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)> {
        match self.root.intersect(ray, t_min, t_max) {
            Some((t, obj_idx)) => Some((t, &self.objects[obj_idx])),
            None => None,
        }
    }
}

impl<T: Bounded<D>, const D: usize> Bvh<T, D> {
    pub fn build(objects: Vec<T>) -> Self {
        let mut sorted: Vec<Vec<usize>> = (0..D).map(|_| (0..objects.len()).collect()).collect();
        sorted.iter_mut().enumerate().for_each(|(i, v)| {
            v.sort_by(|a, b| {
                objects[*a].bounds().centroid()[i]
                    .partial_cmp(&objects[*b].bounds().centroid()[i])
                    .unwrap()
            })
        });
        let sorted_slice: Vec<&[_]> = sorted.iter().map(|v| v.as_slice()).collect();

        let root = Self::_build(&objects, &sorted_slice);

        Self { objects, root }
    }

    fn _build(objects: &[T], sorted: &[&[usize]]) -> BvhNode<D> {
        let bounds = sorted[0]
            .iter()
            .map(|i| objects[*i].bounds())
            .reduce(|acc, b| acc.union(&b))
            .expect("No objects to build bounds.");

        match sorted[0].len() {
            0 => panic!("No objects given."),
            1 => BvhNode::Leaf {
                bounds: objects[sorted[0][0]].bounds(),
                data: sorted[0][0],
            },
            2 => BvhNode::Node {
                bounds,
                children: [
                    Box::new(BvhNode::Leaf {
                        bounds: objects[sorted[0][0]].bounds(),
                        data: sorted[0][0],
                    }),
                    Box::new(BvhNode::Leaf {
                        bounds: objects[sorted[0][1]].bounds(),
                        data: sorted[0][1],
                    }),
                ],
            },
            _ => {
                const NUM_BUCKETS: usize = 16;
                let mut scores = [[0.0; NUM_BUCKETS]; D];
                let sa_bounds = bounds.surface_area();

                for axis in 0..D {
                    for bucket in 0..NUM_BUCKETS {
                        let cut_location = bounds.min[axis]
                            + bounds.length(axis)
                                * ((bucket + 1) as f32 / (NUM_BUCKETS + 1) as f32);
                        let (left, right) =
                            sorted[axis].split_at(sorted[axis].partition_point(|o| {
                                objects[*o].bounds().centroid()[axis] < cut_location
                            }));

                        // Bad cut location one of the sides has no shapes.
                        if left.len() == 0 || right.len() == 0 {
                            scores[axis][bucket] = f32::INFINITY;
                            continue;
                        }

                        let left_sa = left
                            .iter()
                            .map(|o| objects[*o].bounds())
                            .reduce(|acc, b| acc.union(&b))
                            .expect("No left bounds!")
                            .surface_area();
                        let right_sa = right
                            .iter()
                            .map(|o| objects[*o].bounds())
                            .reduce(|acc, b| acc.union(&b))
                            .expect("No right bounds!")
                            .surface_area();

                        let left_ratio = left_sa / sa_bounds;
                        let right_ratio = right_sa / sa_bounds;

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

                // Split the shapes into two groups
                let split_index = if scores
                    .iter()
                    .all(|dim| dim.iter().all(|v| v == &f32::INFINITY))
                {
                    // Bound centroids are all on top of eachother. Split in the middle.
                    sorted[min_axis].len() / 2
                } else {
                    let cut_location = bounds.min[min_axis]
                        + bounds.length(min_axis)
                            * ((min_bucket + 1) as f32 / (NUM_BUCKETS + 1) as f32);

                    sorted[min_axis].partition_point(|o| {
                        objects[*o].bounds().centroid()[min_axis] < cut_location
                    })
                };

                let mut left: Vec<Vec<usize>> = (0..D).map(|_| sorted[min_axis].to_vec()).collect();
                let mut right: Vec<Vec<usize>> =
                    left.iter_mut().map(|v| v.split_off(split_index)).collect();

                // Re-sort the axes
                left.iter_mut().enumerate().for_each(|(i, v)| {
                    v.sort_by(|a, b| {
                        objects[*a].bounds().centroid()[i]
                            .partial_cmp(&objects[*b].bounds().centroid()[i])
                            .unwrap()
                    })
                });
                right.iter_mut().enumerate().for_each(|(i, v)| {
                    v.sort_by(|a, b| {
                        objects[*a].bounds().centroid()[i]
                            .partial_cmp(&objects[*b].bounds().centroid()[i])
                            .unwrap()
                    })
                });

                let left_slice: Vec<&[_]> = left.iter().map(|v| v.as_slice()).collect();
                let right_slice: Vec<&[_]> = right.iter().map(|v| v.as_slice()).collect();

                let left_node = Self::_build(objects, &left_slice);
                let right_node = Self::_build(objects, &right_slice);

                BvhNode::Node {
                    bounds,
                    children: [Box::new(left_node), Box::new(right_node)],
                }
            }
        }
    }
}

impl<'a, T: Bounded<D>, const D: usize> IntoIterator for &'a Bvh<T, D> {
    type Item = BvhIteratorNode<'a, T, D>;
    type IntoIter = BvhIterator<'a, T, D>;

    fn into_iter(self) -> Self::IntoIter {
        BvhIterator {
            bvh: self,
            stack: vec![&self.root],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BvhIteratorNode<'a, T: Bounded<D>, const D: usize> {
    Node(&'a Bounds<D>),
    Leaf(&'a T),
}

#[derive(Debug, Clone)]
pub struct BvhIterator<'a, T: Bounded<D>, const D: usize> {
    bvh: &'a Bvh<T, D>,
    stack: Vec<&'a BvhNode<D>>,
}

impl<'a, T: Bounded<D>, const D: usize> Iterator for BvhIterator<'a, T, D> {
    type Item = BvhIteratorNode<'a, T, D>;

    fn next(&mut self) -> Option<Self::Item> {
        let node = match self.stack.pop() {
            Some(node) => node,
            None => return None,
        };

        match node {
            BvhNode::Leaf { data, .. } => {
                let object = &self.bvh.objects[*data];
                Some(BvhIteratorNode::Leaf(object))
            }
            BvhNode::Node { bounds, children } => {
                self.stack.push(&children[1]);
                self.stack.push(&children[0]);

                Some(BvhIteratorNode::Node(&bounds))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn bounds_length_2d() {
        let b = Bounds::new([0.0, 0.0], [1.2, 1.3]);
        assert_eq!(b.length(0), 1.2);
        assert_eq!(b.length(1), 1.3);
    }

    #[test]
    fn bounds_length_3d() {
        let b = Bounds::new([0.0, 0.0, 0.0], [1.2, 1.3, 1.4]);
        assert_eq!(b.length(0), 1.2);
        assert_eq!(b.length(1), 1.3);
        assert_eq!(b.length(2), 1.4);
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
    fn bounds_intersect_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let r1 = ([0.5, 2.0], [0.0, -1.0]); // intersects
        let r2 = ([0.5, 0.5], [0.0, 1.0]); // intersects from inside
        let r3 = ([0.5, 2.0], [1.0, 0.0]); // miss
        let r4 = ([0.5, 2.0], [0.0, 1.0]); // miss
        assert_eq!(a.intersect(&r1, f32::MIN, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.intersect(&r2, f32::MIN, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.intersect(&r3, f32::MIN, f32::MAX), None);
        assert_eq!(a.intersect(&r4, f32::MIN, f32::MAX), None);

        let b = Bounds::new([3.0, -2.0], [4.0, -1.0]);
        let r = ([4.5, -2.5], [1.0, 1.0]);
        assert_eq!(b.intersect(&r, f32::MIN, f32::MAX), None);

        let b = Bounds::new([4.5, 0.0], [5.5, 1.0]);
        let r = ([4.5, -2.5], [1.0, 1.0]);
        assert_eq!(b.intersect(&r, f32::MIN, f32::MAX), None);
    }

    #[test]
    fn bounds_intersect_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let r1 = ([0.5, 2.0, 0.5], [0.0, -1.0, 0.0]); // intersects
        let r2 = ([0.5, 0.5, 0.5], [0.0, 1.0, 0.0]); // intersects from inside
        let r3 = ([0.5, 2.0, 0.5], [1.0, 0.0, 0.0]); // miss
        let r4 = ([0.5, 2.0, 0.5], [0.0, 1.0, 0.0]); // miss
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

        let r1 = ([3.5, -2.5], [0.0, 1.0]); // hits 0
        let intersection = bvh.intersect(&r1, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 0.5);

        let r2 = ([6.5, 0.5], [-1.0, 0.0]); // hits 1
        let intersection = bvh.intersect(&r2, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 1.0);

        let r3 = ([4.5, -2.5], [1.0, 1.0]); // Hits bounding box but not a leaf
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

        let r1 = ([3.5, 5.0, 3.5], [0.0, -1.0, 0.0]); // hits 0
        let intersection = bvh.intersect(&r1, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 1.0);

        let r2 = ([0.5, 0.5, -0.5], [0.0, 0.0, 1.0]); // hits 1
        let intersection = bvh.intersect(&r2, f32::MIN, f32::MAX);
        assert!(intersection.is_some());
        assert_eq!(intersection.unwrap().0, 0.5);

        let r3 = ([0.5, -0.5, 3.5], [0.0, 1.0, 0.0]); // Hits bounding box but not a leaf
        let intersection = bvh.intersect(&r3, f32::MIN, f32::MAX);
        assert!(intersection.is_none());
    }

    #[test]
    fn bvh_iter_2d() {
        let objects = vec![
            Bounds::new([3.0, -2.0], [4.0, -1.0]), // 0
            Bounds::new([4.5, 0.0], [5.5, 1.0]),   // 1
            Bounds::new([0.0, 1.0], [1.0, 2.0]),   // 2
            Bounds::new([1.5, -0.5], [2.5, 0.5]),  // 3
        ];
        let bvh = Bvh2d::build(objects);

        let mut iter = bvh.into_iter();
        match iter.next().unwrap() {
            BvhIteratorNode::Node(node) => assert_eq!(node.bounds(), bvh.bounds()),
            _ => panic!("Expected node"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Node(node) => {
                assert_eq!(node.bounds(), Bounds::new([0.0, -0.5], [2.5, 2.0]))
            }
            _ => panic!("Expected node"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Leaf(leaf) => {
                assert_eq!(leaf.bounds(), Bounds::new([0.0, 1.0], [1.0, 2.0]))
            }
            _ => panic!("Expected leaf"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Leaf(leaf) => {
                assert_eq!(leaf.bounds(), Bounds::new([1.5, -0.5], [2.5, 0.5]))
            }
            _ => panic!("Expected leaf"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Node(node) => {
                assert_eq!(node.bounds(), Bounds::new([3.0, -2.0], [5.5, 1.0]))
            }
            _ => panic!("Expected node"),
        }
    }

    #[test]
    fn bvh_iter_3d() {
        let objects = vec![
            Bounds::new([3.0, 3.0, 3.0], [4.0, 4.0, 4.0]), // 0
            Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), // 1
            Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]), // 2
            Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]), // 3
        ];
        let bvh = Bvh3d::build(objects);

        let mut iter = bvh.into_iter();
        match iter.next().unwrap() {
            BvhIteratorNode::Node(node) => assert_eq!(node.bounds(), bvh.bounds()),
            _ => panic!("Expected node"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Node(node) => {
                assert_eq!(node.bounds(), Bounds::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]))
            }
            _ => panic!("Expected node"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Leaf(leaf) => {
                assert_eq!(leaf.bounds(), Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]))
            }
            _ => panic!("Expected leaf"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Leaf(leaf) => {
                assert_eq!(leaf.bounds(), Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]))
            }
            _ => panic!("Expected leaf"),
        }
        match iter.next().unwrap() {
            BvhIteratorNode::Node(node) => {
                assert_eq!(node.bounds(), Bounds::new([2.0, 2.0, 2.0], [4.0, 4.0, 4.0]))
            }
            _ => panic!("Expected node"),
        }
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
        let r = ([-0.5, 49.5], [1.0, 0.0]);
        b.iter(|| bvh.intersect(&r, f32::MIN, f32::MAX));
    }

    #[bench]
    fn bench_bvh_intersect_3d(b: &mut Bencher) {
        let objects: Vec<Bounds<3>> = (0..100)
            .map(|i| Bounds::new([i as f32; 3], [(i + 1) as f32; 3]))
            .collect();
        let bvh = Bvh3d::build(objects);
        let r = ([-0.5, 49.5, 49.5], [1.0, 0.0, 0.0]);
        b.iter(|| bvh.intersect(&r, f32::MIN, f32::MAX));
    }

    #[bench]
    fn bench_bvh_iter_2d(b: &mut Bencher) {
        let objects: Vec<Bounds<2>> = (0..100)
            .map(|i| Bounds::new([i as f32; 2], [(i + 1) as f32; 2]))
            .collect();
        let bvh = Bvh2d::build(objects);
        b.iter(|| bvh.into_iter().map(|_| ()).count());
    }

    #[bench]
    fn bench_bvh_iter_3d(b: &mut Bencher) {
        let objects: Vec<Bounds<3>> = (0..100)
            .map(|i| Bounds::new([i as f32; 3], [(i + 1) as f32; 3]))
            .collect();
        let bvh = Bvh3d::build(objects);
        b.iter(|| bvh.into_iter().map(|_| ()).count());
    }
}
