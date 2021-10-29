#[allow(unused_imports)]
use crate::bvh::Bvh;
use std::ops::{Index, IndexMut};
/// Trait definition for a DIM-dimensional Axis-Aligned Bounding Box.
pub trait BoundingBox:
    Copy + RayHittable<Self> + BoundsHittable<Self> + PointHittable<Self>
{
    const DIM: usize;
    type Vector: Index<usize, Output = f32> + IndexMut<usize, Output = f32>;
    type Ray;

    fn min(&self) -> &Self::Vector;
    fn max(&self) -> &Self::Vector;
    fn min_mut(&mut self) -> &mut Self::Vector;
    fn max_mut(&mut self) -> &mut Self::Vector;

    fn shape(&self) -> Self::Vector;
    fn axis_length(&self, axis: usize) -> f32;
    fn centroid(&self) -> Self::Vector;
    fn surface_area(&self) -> f32;
    fn union(&self, other: &Self) -> Self;
}

/// Trait definition for an object that can be contained in a [BoundingBox].
///
/// This trait is required to be implemented in order to use Bvh.
pub trait Bounded<B: BoundingBox> {
    fn bounds(&self) -> B;
}

/// Trait definition for an object that can be intersected by a [Ray](BoundingBox::Ray).
///
/// This trait is NOT required to be implemented in order to build a Bvh.
/// It is reccommended that you implement this trait as it gives you
/// access to [Bvh::query_ray_exact] and [Bvh::ray_hit].
pub trait RayHittable<B: BoundingBox>: Bounded<B> {
    type Item: Copy;

    fn ray_hit(&self, ray: &B::Ray, t_min: f32, t_max: f32) -> Option<(f32, Self::Item)>;
}

/// Trait definition for an object that can be intersected by a [BoundingBox].
///
/// This trait is NOT required to be implemented in order to build a Bvh.
/// It is reccommended that you implement this trait as it gives you
/// access to [Bvh::query_bounds_exact].
pub trait BoundsHittable<B: BoundingBox>: Bounded<B> {
    fn bounds_hit(&self, bounds: &B) -> bool;
}

/// Trait definition for an object that can be intersected by a [Point](BoundingBox::Vector)
///
/// This trait is NOT required to be implemented in order to build a Bvh.
/// It is reccommended that you implement this trait as it gives you
/// access to [Bvh::query_point_exact].
pub trait PointHittable<B: BoundingBox>: Bounded<B> {
    fn point_hit(&self, point: &B::Vector) -> bool;
}
