use crate::bounds::{Bounds2, Bounds3, Bounds3A};
use crate::{Ray2, Ray3, Ray3A};

use glam::{Vec2, Vec3, Vec3A};

pub trait Bounded2: Clone {
    const DIM: usize = 2;
    fn bounds(&self) -> Bounds2;
}
pub trait Ray2Hit: Bounded2 {
    type Item: Bounded2;
    fn ray_hit(&self, ray: &Ray2, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)>;
}
pub trait Bounds2Hit: Bounded2 {
    fn bounds_hit(&self, bounds: &Bounds2) -> bool;
}
pub trait Point2Hit: Bounded2 {
    fn point_hit(&self, point: &Vec2) -> bool;
}

pub trait Bounded3: Clone {
    const DIM: usize = 3;
    fn bounds(&self) -> Bounds3;
}
pub trait Ray3Hit: Bounded3 {
    type Item: Bounded3;
    fn ray_hit(&self, ray: &Ray3, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)>;
}
pub trait Bounds3Hit: Bounded3 {
    fn bounds_hit(&self, bounds: &Bounds3) -> bool;
}
pub trait Point3Hit: Bounded3 {
    fn point_hit(&self, point: &Vec3) -> bool;
}

pub trait Bounded3A: Clone {
    const DIM: usize = 3;
    fn bounds(&self) -> Bounds3A;
}
pub trait Ray3AHit: Bounded3A {
    type Item: Bounded3A;
    fn ray_hit(&self, ray: &Ray3A, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)>;
}
pub trait Bounds3AHit: Bounded3A {
    fn bounds_hit(&self, bounds: &Bounds3A) -> bool;
}
pub trait Point3AHit: Bounded3A {
    fn point_hit(&self, point: &Vec3A) -> bool;
}
