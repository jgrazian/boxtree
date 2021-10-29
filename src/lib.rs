pub mod bounds;
pub mod bvh;
mod iter;
mod split;
pub mod traits;

pub use bounds::*;
pub use bvh::*;
pub use ray::*;
pub use split::*;
pub use traits::*;

pub extern crate glam;

macro_rules! impl_ray {
    ($name:ident, $vec:ty, $dim:expr, $comment:expr) => {
        #[doc=$comment]
        #[derive(Clone, Copy, Debug)]
        pub struct $name {
            pub origin: $vec,
            pub direction: $vec,
        }

        impl $name {
            /// Creates a new ray at `origin` moving `direction`.
            pub fn new(origin: $vec, direction: $vec) -> Self {
                Self { origin, direction }
            }

            /// Where the ray is in space a given time `t`.
            ///
            /// Calculated as `origin + direction*t`.
            pub fn at(&self, t: f32) -> $vec {
                self.origin + self.direction * t
            }
        }

        impl From<($vec, $vec)> for $name {
            fn from(tuple: ($vec, $vec)) -> Self {
                Self::new(tuple.0, tuple.1)
            }
        }

        impl From<([f32; $dim], [f32; $dim])> for $name {
            fn from(tuple: ([f32; $dim], [f32; $dim])) -> Self {
                Self::new(tuple.0.into(), tuple.1.into())
            }
        }
    };
}

pub mod ray {
    impl_ray!(
        Ray2,
        glam::Vec2,
        2,
        "A 2-dimensional ray using [Vec2](glam::Vec2)."
    );
    impl_ray!(
        Ray3,
        glam::Vec3,
        3,
        "A 3-dimensional ray using [Vec3](glam::Vec3)."
    );
    impl_ray!(
        Ray3A,
        glam::Vec3A,
        3,
        "A 3-dimensional ray using [Vec3A](glam::Vec3A)."
    );
}
