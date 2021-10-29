use std::ops::Index;

use crate::bounds::{Aabb2, Aabb3, Aabb3A};
use crate::iter::BvhLeafIterator;
use crate::split::SplitMethod;
use crate::traits::*;

pub type Bvh2<T> = Bvh<Aabb2, T, 2>;
pub type Bvh3<T> = Bvh<Aabb3, T, 2>;
pub type Bvh3A<T> = Bvh<Aabb3A, T, 2>;

/// Can be used to index into a Bvh and retrieve a reference to an object there.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BvhObjKey(usize);

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct BvhNodeKey(pub(crate) usize);

/// A Bvh Node.
/// A Node will always have 2 children.
/// Leaves contain up to `N` children.
#[derive(Clone, Debug)]
pub(crate) enum BvhNode<B: BoundingBox, const N: usize> {
    Node {
        bounds: B,
        parent: Option<BvhNodeKey>,
        children: [BvhNodeKey; 2],
    },
    Leaf {
        bounds: B,
        parent: Option<BvhNodeKey>,
        objects: [Option<BvhObjKey>; N],
    },
}

/// Bounding Volume Hierarchy.
///
/// `N` is the maximum number of objects that a single leaf node can hold.
///
/// Internally Bvh holds all of its [Bounded] objects in a [Vec]
/// and all nodes that make up the tree in another [Vec]. Each node holds a reference to its
/// parent node and either its child nodes or child objects if it is a leaf.
#[derive(Clone, Debug)]
pub struct Bvh<B: BoundingBox, T: Bounded<B>, const N: usize> {
    objects: Vec<T>,
    pub(crate) nodes: Vec<BvhNode<B, N>>,
}

impl<'a, B: BoundingBox, T: Bounded<B> + Clone, const N: usize> Bvh<B, T, N> {
    /// Constructs a Bvh from a slice of `Bounded + Clone` objects.
    ///
    /// The build algorithm is a stack based iteration.
    pub fn build(objects: &[T], split_method: SplitMethod) -> Self {
        let mut out_objects: Vec<T> = Vec::with_capacity(objects.len());
        let centroids = objects
            .iter()
            .map(|obj| obj.bounds().centroid())
            .collect::<Vec<_>>();
        let mut indexes = (0..objects.len()).collect::<Vec<_>>();

        let mut nodes = Vec::with_capacity(Self::STACK_SIZE);
        let mut stack: Vec<(Option<BvhNodeKey>, usize, &mut [usize])> =
            vec![(None, 0, &mut indexes)];

        while let Some((parent, child_pos, idxs)) = stack.pop() {
            match idxs.len() {
                0 => (),
                n if n <= N => {
                    let bounds = Self::grow_bounds(objects, idxs);

                    let mut leaf_objects = [None; N];
                    for (i, idx) in idxs.iter().enumerate() {
                        out_objects.push(objects[*idx].clone());
                        leaf_objects[i] = Some(BvhObjKey(out_objects.len() - 1));
                    }

                    nodes.push(BvhNode::Leaf {
                        bounds,
                        parent,
                        objects: leaf_objects,
                    });
                    let node_id = BvhNodeKey(nodes.len() - 1);

                    if let Some(parent_key) = parent {
                        match nodes[parent_key.0] {
                            BvhNode::Node {
                                ref mut children, ..
                            } => children[child_pos] = node_id,
                            _ => panic!("Expected node found leaf."),
                        }
                    }
                }
                _ => {
                    let bounds = Self::grow_bounds(objects, idxs);

                    nodes.push(BvhNode::Node {
                        bounds,
                        parent,
                        children: [BvhNodeKey::default(); 2],
                    });
                    let node_id = BvhNodeKey(nodes.len() - 1);

                    if let Some(parent_key) = parent {
                        match nodes[parent_key.0] {
                            BvhNode::Node {
                                ref mut children, ..
                            } => children[child_pos] = node_id,
                            _ => panic!("Expected node found leaf."),
                        }
                    }

                    let (left, right) = split_method.split(&bounds, &objects, &centroids, idxs);
                    stack.push((Some(node_id), 1, right));
                    stack.push((Some(node_id), 0, left));
                }
            }
        }

        Self {
            objects: out_objects,
            nodes,
        }
    }

    fn grow_bounds(objects: &[T], indexes: &[usize]) -> B {
        indexes
            .iter()
            .map(|i| objects[*i].bounds())
            .reduce(|acc, b| acc.union(&b))
            .expect("No objects to build bounds.")
    }
}

impl<B: BoundingBox, T: Bounded<B>, const N: usize> Bounded<B> for Bvh<B, T, N> {
    fn bounds(&self) -> B {
        match self.nodes[0] {
            BvhNode::Node { bounds, .. } => bounds,
            BvhNode::Leaf { bounds, .. } => bounds,
        }
    }
}

impl<'a, B: BoundingBox, T: Bounded<B>, const N: usize> Bvh<B, T, N> {
    const STACK_SIZE: usize = 32;

    /// Creates an iterator over the objects contained within the Bvh in normal DFS traversal order.
    /// If `predicate` is `false` that node and all children will be skipped during iteration.
    pub fn iter_objects_pred<F: 'a + Fn(&B) -> bool>(
        &'a self,
        predicate: F,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        BvhLeafIterator::new(&self, predicate)
    }
    /// Iterates over all objects in Bvh in DFS order.
    pub fn iter_objects(&'a self) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        self.iter_objects_pred(|_| true)
    }
    /// Creates an iterator of references to all objects who's parent leaf's
    /// bounding box intersects with the query ray between times `t_min` and `t_max`.
    ///
    /// See [query_ray_exact](Self::query_ray_exact) for intersections
    /// between the query ray and objects directly.
    pub fn query_ray(
        &'a self,
        ray: &'a B::Ray,
        t_min: f32,
        t_max: f32,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let predicate = move |bounds: &B| bounds.ray_hit(&ray, t_min, t_max).is_some();

        self.iter_objects_pred(predicate)
    }
    /// Creates an iterator of references to all objects who's parent leaf's
    /// bounding box intersects with the query bounding box.
    ///
    /// See [query_bounds_exact](Self::query_bounds_exact) for intersection
    /// between bounds and objects directly.
    pub fn query_bounds(&'a self, bounds: &'a B) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let predicate = move |bounds2: &B| bounds.bounds_hit(&bounds2);

        self.iter_objects_pred(predicate)
    }
    /// Creates an iterator of references to all objects who's parent leaf's
    /// bounding box intersects with the query bounding box.
    ///
    /// See [query_point_exact](Self::query_point_exact) for intersection
    /// between point and objects directly.
    pub fn query_point(
        &'a self,
        point: &'a B::Vector,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        let predicate = move |bounds: &B| bounds.point_hit(&point);

        self.iter_objects_pred(predicate)
    }
}

impl<'a, B: BoundingBox, T: RayHittable<B>, const N: usize> Bvh<B, T, N> {
    /// Creates an iterator of references to all objects who intersects with the query ray.
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
    /// Creates an iterator of references to all objects who intersects with the query bounding box.
    pub fn query_bounds_exact(
        &'a self,
        bounds: &'a B,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        self.query_bounds(bounds)
            .filter(move |(_, obj)| obj.bounds_hit(bounds))
    }
}

impl<'a, B: BoundingBox, T: PointHittable<B>, const N: usize> Bvh<B, T, N> {
    /// Creates an iterator of references to all objects who intersects with the query point.
    pub fn query_point_exact(
        &'a self,
        point: &'a B::Vector,
    ) -> impl Iterator<Item = (BvhObjKey, &T)> + '_ {
        self.query_point(point)
            .filter(move |(_, obj)| obj.point_hit(point))
    }
}

impl<B: BoundingBox, T: RayHittable<B>, const N: usize> RayHittable<B> for Bvh<B, T, N> {
    type Item = T::Item;

    /// Find the closest object hit by a ray.
    ///
    /// You could use [query_ray_exact](Self::query_ray_exact) to perform this function as well
    /// but this method is optimized better.
    fn ray_hit(&self, ray: &B::Ray, t_min: f32, t_max: f32) -> Option<(f32, Self::Item)> {
        let mut stack = Vec::with_capacity(Self::STACK_SIZE);
        stack.push(&self.nodes[0]);
        let mut result = None;

        while let Some(node) = stack.pop() {
            let t_max = result.map_or(t_max, |(t, _)| t);
            match node {
                BvhNode::Node {
                    bounds, children, ..
                } => {
                    if bounds.ray_hit(ray, t_min, t_max).is_some() {
                        for child_key in children.iter().rev() {
                            stack.push(&self.nodes[child_key.0]);
                        }
                    }
                }
                BvhNode::Leaf {
                    bounds, objects, ..
                } => {
                    if bounds.ray_hit(ray, t_min, t_max).is_some() {
                        for obj_key in objects.iter() {
                            if let Some(key) = obj_key {
                                let obj = &self.objects[key.0];
                                if let Some((t, item)) = obj.ray_hit(ray, t_min, t_max) {
                                    result = Some((t, item));
                                }
                            }
                        }
                    }
                }
            }
        }
        result
    }
}

impl<B: BoundingBox, T: Bounded<B>, const N: usize> Index<BvhObjKey> for Bvh<B, T, N> {
    type Output = T;

    fn index(&self, index: BvhObjKey) -> &Self::Output {
        &self.objects[index.0]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// │╔══╤═╗ ╔══╤═╗
    /// │╟─┐└─╢ ╟─┐└─╢
    /// │╚═╧══╝ ╚═╧══╝
    /// ┼───────────────
    fn setup_2d() -> Bvh2<Aabb2> {
        let objects = vec![
            ([10.0, 0.0], [11.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([11.0, 1.0], [12.0, 2.0]).into(),
            ([1.0, 1.0], [2.0, 2.0]).into(),
        ];
        Bvh2::build(&objects, SplitMethod::SurfaceAreaHeuristic)
    }

    fn setup_3d() -> Bvh3<Aabb3> {
        let objects = vec![
            ([10.0, 0.0, 0.0], [11.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([11.0, 1.0, 1.0], [12.0, 2.0, 2.0]).into(),
            ([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]).into(),
        ];
        Bvh3::build(&objects, SplitMethod::SurfaceAreaHeuristic)
    }

    #[test]
    fn build_2d() {
        let bvh = setup_2d();
        assert_eq!(bvh.bounds(), ([0.0, 0.0], [12.0, 2.0]).into());

        let objects: Vec<Aabb2> = vec![
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
            ([0.0, 0.0], [1.0, 1.0]).into(),
        ];
        let bvh = Bvh2::build(&objects, SplitMethod::SurfaceAreaHeuristic);
        assert_eq!(bvh.bounds(), ([0.0, 0.0], [1.0, 1.0]).into());
    }

    #[test]
    fn build_3d() {
        let bvh = setup_3d();
        assert_eq!(bvh.bounds(), ([0.0, 0.0, 0.0], [12.0, 2.0, 2.0]).into());

        let objects: Vec<Aabb3> = vec![
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into(),
        ];
        let bvh = Bvh3::build(&objects, SplitMethod::SurfaceAreaHeuristic);
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
            (10.0, BvhObjKey(2), &bvh.objects[2])
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
            (10.0, BvhObjKey(2), &bvh.objects[2])
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

        let b = ([0.0, 0.0], [11.5, 0.5]).into();
        let mut bounds_hits = bvh.query_bounds_exact(&b);
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(0), &bvh.objects[0]));
        assert_eq!(bounds_hits.next().unwrap(), (BvhObjKey(2), &bvh.objects[2]));
        assert_eq!(bounds_hits.next(), None);
    }

    #[test]
    fn query_bounds_exact_3d() {
        let bvh = setup_3d();

        let b = ([0.0, 0.0, 0.0], [11.5, 0.5, 0.5]).into();
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

    #[test]
    fn ray_hit_2d() {
        let bvh = setup_2d();

        let r = ([-1.0, 0.5], [1.0, 0.0]).into();
        let hit = bvh.ray_hit(&r, 0.0, f32::MAX);
        assert_eq!(hit, Some((1.0, bvh.objects[0])));
    }

    #[test]
    fn ray_hit_3d() {
        let bvh = setup_3d();

        let r = ([-1.0, 0.5, 0.5], [1.0, 0.0, 0.0]).into();
        let hit = bvh.ray_hit(&r, 0.0, f32::MAX);
        assert_eq!(hit, Some((1.0, bvh.objects[0])));
    }
}
