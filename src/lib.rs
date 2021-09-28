mod bounds;
mod bvh;
mod iter;
mod traits;

pub use bounds::{Bounds2, Bounds3, Bounds3A};
pub use bvh::{Bvh2, Bvh3, Bvh3A};
pub use traits::*;

pub use glam::{Vec2, Vec3, Vec3A};

pub type Ray2 = Ray<Vec2>;
pub type Ray3 = Ray<Vec3>;
pub type Ray3A = Ray<Vec3A>;

/// A ray in D dimensional space.
/// Starts at `origin` and goes in the direction `direction`.
#[derive(Clone, Copy, Debug)]
pub struct Ray<V> {
    origin: V,
    direction: V,
}

impl<V> Ray<V> {
    fn new(origin: V, direction: V) -> Self {
        Self { origin, direction }
    }
}

impl<V> From<(V, V)> for Ray<V> {
    fn from(tuple: (V, V)) -> Self {
        Self::new(tuple.0, tuple.1)
    }
}

// Vec2
impl Ray2 {
    #[allow(dead_code)]
    fn at(&self, t: f32) -> Vec2 {
        self.origin + self.direction * t
    }
}
impl From<([f32; 2], [f32; 2])> for Ray2 {
    fn from(tuple: ([f32; 2], [f32; 2])) -> Self {
        Self::new(tuple.0.into(), tuple.1.into())
    }
}

// Vec3
impl Ray3 {
    #[allow(dead_code)]
    fn at(&self, t: f32) -> Vec3 {
        self.origin + self.direction * t
    }
}
impl From<([f32; 3], [f32; 3])> for Ray3 {
    fn from(tuple: ([f32; 3], [f32; 3])) -> Self {
        Self::new(tuple.0.into(), tuple.1.into())
    }
}

// Vec3A
impl Ray3A {
    #[allow(dead_code)]
    fn at(&self, t: f32) -> Vec3A {
        self.origin + self.direction * t
    }
}
impl From<([f32; 3], [f32; 3])> for Ray3A {
    fn from(tuple: ([f32; 3], [f32; 3])) -> Self {
        Self::new(tuple.0.into(), tuple.1.into())
    }
}
