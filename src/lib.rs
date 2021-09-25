mod bounds;
mod bvh;
mod iter;
mod traits;

pub use bounds::*;
pub use bvh::{Bvh2, Bvh3, Bvh3A};
pub use traits::*;

pub use glam::{Vec2, Vec3, Vec3A};

#[derive(Clone, Copy, Debug)]
pub struct Ray2 {
    origin: Vec2,
    direction: Vec2,
}

#[derive(Clone, Copy, Debug)]
pub struct Ray3 {
    origin: Vec3,
    direction: Vec3,
}

#[derive(Clone, Copy, Debug)]
pub struct Ray3A {
    origin: Vec3A,
    direction: Vec3A,
}

macro_rules! impl_ray {
    ($dim:tt, $ray:ty, $vec:ty) => {
        impl $ray {
            fn new(origin: $vec, direction: $vec) -> Self {
                Self { origin, direction }
            }

            pub fn at(&self, t: f32) -> $vec {
                self.origin + self.direction * t
            }
        }
        impl From<($vec, $vec)> for $ray {
            fn from(tuple: ($vec, $vec)) -> Self {
                Self::new(tuple.0, tuple.1)
            }
        }
        impl From<([f32; $dim], [f32; $dim])> for $ray {
            fn from(tuple: ([f32; $dim], [f32; $dim])) -> Self {
                Self::new(tuple.0.into(), tuple.1.into())
            }
        }
    };
}

impl_ray!(2, Ray2, Vec2);
impl_ray!(3, Ray3, Vec3);
impl_ray!(3, Ray3A, Vec3A);
