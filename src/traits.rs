use crate::bounds::Bounds;
use crate::{Ray, Vector};

pub trait Bounded<const D: usize>: Clone {
    type Bound: Bounded<D>;

    fn bounds(&self) -> Bounds<D>;
}

pub trait RayHittable<const D: usize>: Bounded<D> {
    fn ray_hit(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)>;
}

pub trait BoundsHittable<const D: usize>: Bounded<D> {
    fn bounds_hit(&self, bounds: &Bounds<D>) -> bool;
}

pub trait PointHittable<const D: usize>: Bounded<D> {
    fn point_hit(&self, point: &Vector<D>) -> bool;
}
