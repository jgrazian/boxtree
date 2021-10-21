use crate::traits::*;
use crate::{Ray2, Ray3, Ray3A};

use glam::{Vec2, Vec3, Vec3A};

pub type Aabb2 = Aabb<Vec2>;
pub type Aabb3 = Aabb<Vec3>;
pub type Aabb3A = Aabb<Vec3A>;

/// A bounding box in D dimensional space.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Aabb<V> {
    min: V,
    max: V,
}

macro_rules! impl_bounds {
    ($dim:expr, $vec:ident, $ray:ident) => {
        impl Aabb<$vec> {
            pub fn new(min: $vec, max: $vec) -> Self {
                Aabb { min, max }
            }
        }

        impl From<($vec, $vec)> for Aabb<$vec> {
            fn from(tuple: ($vec, $vec)) -> Self {
                Self::new(tuple.0, tuple.1)
            }
        }
        impl From<([f32; $dim], [f32; $dim])> for Aabb<$vec> {
            fn from(tuple: ([f32; $dim], [f32; $dim])) -> Self {
                Self::new(tuple.0.into(), tuple.1.into())
            }
        }

        impl BoundingBox for Aabb<$vec> {
            const DIM: usize = $dim;
            type Vector = $vec;
            type Ray = $ray;

            fn min(&self) -> &Self::Vector {
                &self.min
            }
            fn max(&self) -> &Self::Vector {
                &self.max
            }
            fn min_mut(&mut self) -> &mut Self::Vector {
                &mut self.min
            }
            fn max_mut(&mut self) -> &mut Self::Vector {
                &mut self.max
            }

            fn shape(&self) -> Self::Vector {
                self.max - self.min
            }

            fn axis_length(&self, axis: usize) -> f32 {
                self.max[axis] - self.min[axis]
            }

            fn centroid(&self) -> Self::Vector {
                (self.min + self.max) * 0.5
            }

            fn surface_area(&self) -> f32 {
                let shape = self.shape();
                let mut area = 0.0;
                for i in 0..Self::DIM {
                    for j in (i + 1)..Self::DIM {
                        area += shape[i] * shape[j];
                    }
                }
                area * (Self::DIM - 1) as f32
            }

            fn union(&self, other: &Self) -> Self {
                Self::new(self.min.min(other.min), self.max.max(other.max))
            }
        }

        impl Bounded<Aabb<$vec>> for Aabb<$vec> {
            fn bounds(&self) -> Self {
                *self
            }
        }

        impl RayHittable<Aabb<$vec>> for Aabb<$vec> {
            type Item = Self;

            fn ray_hit(
                &self,
                ray: &<Self as BoundingBox>::Ray,
                t_min: f32,
                t_max: f32,
            ) -> Option<(f32, Self::Item)> {
                let inv_d = ray.direction.recip();
                let t0__ = (self.min - ray.origin) * inv_d;
                let t1__ = (self.max - ray.origin) * inv_d;

                let mask = inv_d.cmplt(<Self as BoundingBox>::Vector::ZERO);
                let t0_ = <Self as BoundingBox>::Vector::select(mask, t1__, t0__);
                let t1_ = <Self as BoundingBox>::Vector::select(mask, t0__, t1__);

                let t0 = t0_.max_element().max(t_min);
                let t1 = t1_.min_element().min(t_max);

                if t1 <= t0 {
                    return None;
                }

                Some((t0, *self))
            }
        }

        impl BoundsHittable<Aabb<$vec>> for Aabb<$vec> {
            fn bounds_hit(&self, bounds: &Self) -> bool {
                (self.min.cmple(bounds.max) & self.max.cmpge(bounds.min)).all()
            }
        }

        impl PointHittable<Aabb<$vec>> for Aabb<$vec> {
            fn point_hit(&self, point: &<Self as BoundingBox>::Vector) -> bool {
                (self.min.cmple(*point) & self.max.cmpge(*point)).all()
            }
        }
    };
}

impl_bounds!(2, Vec2, Ray2);
impl_bounds!(3, Vec3, Ray3);
impl_bounds!(3, Vec3A, Ray3A);

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn axis_length_2d() {
        let b: Aabb2 = ([0.0, 0.0], [1.2, 1.3]).into();
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
    }

    #[test]
    fn axis_length_3d() {
        let b: Aabb3 = ([0.0, 0.0, 0.0], [1.2, 1.3, 1.4]).into();
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
        assert_eq!(b.axis_length(2), 1.4);
    }

    #[test]
    fn surface_area_2d() {
        let a: Aabb2 = ([0.0, 0.0], [1.0, 1.0]).into();
        assert_eq!(a.surface_area(), 1.0);
    }

    #[test]
    fn surface_area_3d() {
        let a: Aabb3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        assert_eq!(a.surface_area(), 6.0);
    }

    #[test]
    fn union_2d() {
        let a: Aabb2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let b: Aabb2 = ([1.0, 1.0], [2.0, 2.0]).into();
        assert_eq!(a.union(&b), ([0.0, 0.0], [2.0, 2.0]).into());
    }

    #[test]
    fn union_3d() {
        let a: Aabb3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let b: Aabb3 = ([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]).into();
        assert_eq!(a.union(&b), ([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]).into());
    }

    #[test]
    fn ray_hit_2d() {
        let a: Aabb2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let r1 = ([0.5, 2.0], [0.0, -1.0]).into();
        let r2 = ([0.5, 0.5], [0.0, 1.0]).into();
        let r3 = ([0.5, 2.0], [1.0, 0.0]).into();
        let r4 = ([0.5, 2.0], [0.0, 1.0]).into();
        assert_eq!(a.ray_hit(&r1, 0.0, f32::MAX), Some((1.0, a))); // hit
        assert_eq!(a.ray_hit(&r2, 0.0, f32::MAX), Some((0.0, a))); // hit inside
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((-0.5, a)));
        assert_eq!(a.ray_hit(&r3, 0.0, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, 0.0, f32::MAX), None);

        let r = ([0.0, -0.5], [0.0, 1.0]).into();
        assert_eq!(a.ray_hit(&r, 0.0, 0.4), None);
        assert_eq!(a.ray_hit(&r, 1.6, f32::MAX), None);
        assert_eq!(a.ray_hit(&r, 1.4, f32::MAX), Some((1.4, a)));
    }

    #[test]
    fn ray_hit_3d() {
        let a: Aabb3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let r1 = ([0.5, 2.0, 0.5], [0.0, -1.0, 0.0]).into(); // intersects
        let r2 = ([0.5, 0.5, 0.5], [0.0, 1.0, 0.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0, 0.5], [1.0, 0.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0, 0.5], [0.0, 1.0, 0.0]).into(); // miss
        assert_eq!(a.ray_hit(&r1, 0.0, f32::MAX), Some((1.0, a)));
        assert_eq!(a.ray_hit(&r2, 0.0, f32::MAX), Some((0.0, a)));
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((-0.5, a)));
        assert_eq!(a.ray_hit(&r3, 0.0, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, 0.0, f32::MAX), None);

        let r = ([0.5, 0.5, -0.5], [0.0, 0.0, 1.0]).into();
        assert_eq!(a.ray_hit(&r, 0.0, 0.4), None);
        assert_eq!(a.ray_hit(&r, 1.6, f32::MAX), None);
        assert_eq!(a.ray_hit(&r, 1.4, f32::MAX), Some((1.4, a)));
    }

    #[test]
    fn bounds_hit_2d() {
        let a: Aabb2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let b = ([2.0, 2.0], [3.0, 3.0]).into();
        assert!(!a.bounds_hit(&b));
        let c = ([0.5, 0.5], [1.5, 1.5]).into();
        assert!(a.bounds_hit(&c));
    }

    #[test]
    fn bounds_hit_3d() {
        let a: Aabb3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let b = ([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]).into();
        assert!(!a.bounds_hit(&b));
        let c = ([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]).into();
        assert!(a.bounds_hit(&c));
    }

    #[test]
    fn point_hit_2d() {
        let a: Aabb2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let p = [2.0, 2.0].into();
        assert!(!a.point_hit(&p));
        let p = [0.5, 0.5].into();
        assert!(a.point_hit(&p));
    }

    #[test]
    fn point_hit_3d() {
        let a: Aabb3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let p = [2.0, 2.0, 2.0].into();
        assert!(!a.point_hit(&p));
        let p = [0.5, 0.5, 0.5].into();
        assert!(a.point_hit(&p));
    }
}
