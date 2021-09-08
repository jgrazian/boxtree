mod bounds;
mod bvh;
mod iter;
mod traits;

pub use bounds::Bounds;
pub use bvh::{Bvh, Bvh2d, Bvh3d};
pub use traits::*;

/// A D dimensional vector.
type Vector<const D: usize> = [f32; D];

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
