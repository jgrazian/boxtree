use std::ops::{Index, IndexMut};

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

pub trait Bounded<B: BoundingBox> {
    type Item: Bounded<B>;

    fn bounds(&self) -> B;
}

pub trait RayHittable<B: BoundingBox>: Bounded<B> {
    fn ray_hit(&self, ray: &B::Ray, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)>;
}

pub trait BoundsHittable<B: BoundingBox>: Bounded<B> {
    fn bounds_hit(&self, bounds: &B) -> bool;
}

pub trait PointHittable<B: BoundingBox>: Bounded<B> {
    fn point_hit(&self, point: &B::Vector) -> bool;
}
