use crate::traits::*;
use crate::{Ray2, Ray3, Ray3A, Vec2, Vec3, Vec3A};

/// A bounding box in D dimensional space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds2 {
    pub min: Vec2,
    pub max: Vec2,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds3 {
    pub min: Vec3,
    pub max: Vec3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds3A {
    pub min: Vec3A,
    pub max: Vec3A,
}

macro_rules! impl_bounds {
    ($dim:tt, $bounds:ty, $ray:ty, $vec:ty) => {
        impl $bounds {
            const DIM: usize = $dim;

            /// Creates a new $bounds.
            pub fn new(min: $vec, max: $vec) -> Self {
                Self { min, max }
            }

            pub fn shape(&self) -> $vec {
                self.max - self.min
            }

            pub fn axis_length(&self, axis: usize) -> f32 {
                self.max[axis] - self.min[axis]
            }

            pub fn centroid(&self) -> $vec {
                0.5 + (self.min + self.max)
            }

            pub fn surface_area(&self) -> f32 {
                let mut total_area = 0.0;
                for i in 0..Self::DIM {
                    let length = self.axis_length(i);
                    for j in (i + 1)..Self::DIM {
                        let length2 = self.axis_length(j);
                        total_area += length * length2;
                    }
                }
                (Self::DIM - 1) as f32 * total_area
            }

            pub fn union(&self, other: &Self) -> Self {
                Self {
                    min: self.min.min(other.min),
                    max: self.max.max(other.max),
                }
            }
        }

        impl From<($vec, $vec)> for $bounds {
            fn from(tuple: ($vec, $vec)) -> Self {
                Self::new(tuple.0, tuple.1)
            }
        }
        impl From<([f32; $dim], [f32; $dim])> for $bounds {
            fn from(tuple: ([f32; $dim], [f32; $dim])) -> Self {
                Self::new(tuple.0.into(), tuple.1.into())
            }
        }

        impl Bounded<$vec> for $bounds {
            type Item = Self;
            type Bounds = Self;

            fn bounds(&self) -> Self::Bounds {
                *self
            }
        }

        impl RayHittable<$vec, $ray> for $bounds {
            fn ray_hit(&self, ray: &$ray, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)> {
                let mut loc_t_min = t_min;
                let mut loc_t_max = t_max;

                for a in 0..Self::DIM {
                    let inv_d = 1.0 / ray.direction[a];
                    let mut t0 = (self.min[a] - ray.origin[a]) * inv_d;
                    let mut t1 = (self.max[a] - ray.origin[a]) * inv_d;

                    if inv_d < 0.0 {
                        std::mem::swap(&mut t0, &mut t1);
                    }
                    loc_t_min = if t0 > loc_t_min { t0 } else { loc_t_min };
                    loc_t_max = if t1 < loc_t_max { t1 } else { loc_t_max };
                    if loc_t_max <= loc_t_min {
                        return None;
                    }
                }
                match (loc_t_min < 0.0, loc_t_max < 0.0) {
                    (true, true) => None,
                    (true, false) => Some((loc_t_max, &self)),
                    (false, true) => Some((loc_t_min, &self)),
                    (false, false) => Some((loc_t_min, &self)),
                }
            }
        }

        impl BoundsHittable<$vec> for $bounds {
            fn bounds_hit(&self, bounds: &Self::Bounds) -> bool {
                (self.min.cmple(bounds.max) & self.max.cmpge(bounds.min)).all()
            }
        }

        impl PointHittable<$vec> for $bounds {
            fn point_hit(&self, point: &$vec) -> bool {
                (self.min.cmple(*point) & self.max.cmpge(*point)).all()
            }
        }
    };
}

impl_bounds!(2, Bounds2, Ray2, Vec2);
impl_bounds!(3, Bounds3, Ray3, Vec3);
impl_bounds!(3, Bounds3A, Ray3A, Vec3A);

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn axis_length_2d() {
        let b: Bounds2 = ([0.0, 0.0], [1.2, 1.3]).into();
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
    }

    #[test]
    fn axis_length_3d() {
        let b: Bounds3 = ([0.0, 0.0, 0.0], [1.2, 1.3, 1.4]).into();
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
        assert_eq!(b.axis_length(2), 1.4);
    }

    #[test]
    fn surface_area_2d() {
        let a: Bounds2 = ([0.0, 0.0], [1.0, 1.0]).into();
        assert_eq!(a.surface_area(), 1.0);
    }

    #[test]
    fn surface_area_3d() {
        let a: Bounds3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        assert_eq!(a.surface_area(), 6.0);
    }

    #[test]
    fn union_2d() {
        let a: Bounds2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let b: Bounds2 = ([1.0, 1.0], [2.0, 2.0]).into();
        assert_eq!(a.union(&b), ([0.0, 0.0], [2.0, 2.0]).into());
    }

    #[test]
    fn union_3d() {
        let a: Bounds3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let b: Bounds3 = ([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]).into();
        assert_eq!(a.union(&b), ([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]).into());
    }

    #[test]
    fn ray_hit_2d() {
        let a: Bounds2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let r1 = ([0.5, 2.0], [0.0, -1.0]).into(); // intersects
        let r2 = ([0.5, 0.5], [0.0, 1.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0], [1.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0], [0.0, 1.0]).into(); // miss
        assert_eq!(a.ray_hit(&r1, f32::MIN, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.ray_hit(&r3, f32::MIN, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, f32::MIN, f32::MAX), None);

        let b: Bounds2 = ([3.0, -2.0], [4.0, -1.0]).into();
        let r = ([4.5, -2.5], [1.0, 1.0]).into();
        assert_eq!(b.ray_hit(&r, f32::MIN, f32::MAX), None);

        let b: Bounds2 = ([4.5, 0.0], [5.5, 1.0]).into();
        let r = ([4.5, -2.5], [1.0, 1.0]).into();
        assert_eq!(b.ray_hit(&r, f32::MIN, f32::MAX), None);
    }

    #[test]
    fn ray_hit_3d() {
        let a: Bounds3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let r1 = ([0.5, 2.0, 0.5], [0.0, -1.0, 0.0]).into(); // intersects
        let r2 = ([0.5, 0.5, 0.5], [0.0, 1.0, 0.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0, 0.5], [1.0, 0.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0, 0.5], [0.0, 1.0, 0.0]).into(); // miss
        assert_eq!(a.ray_hit(&r1, f32::MIN, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.ray_hit(&r3, f32::MIN, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, f32::MIN, f32::MAX), None);
    }

    #[test]
    fn bounds_hit_2d() {
        let a: Bounds2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let b = ([2.0, 2.0], [3.0, 3.0]).into();
        assert!(!a.bounds_hit(&b));
        let c = ([0.5, 0.5], [1.5, 1.5]).into();
        assert!(a.bounds_hit(&c));
    }

    #[test]
    fn bounds_hit_3d() {
        let a: Bounds3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let b = ([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]).into();
        assert!(!a.bounds_hit(&b));
        let c = ([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]).into();
        assert!(a.bounds_hit(&c));
    }

    #[test]
    fn point_hit_2d() {
        let a: Bounds2 = ([0.0, 0.0], [1.0, 1.0]).into();
        let p = [2.0, 2.0].into();
        assert!(!a.point_hit(&p));
        let p = [0.5, 0.5].into();
        assert!(a.point_hit(&p));
    }

    #[test]
    fn point_hit_3d() {
        let a: Bounds3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let p = [2.0, 2.0, 2.0].into();
        assert!(!a.point_hit(&p));
        let p = [0.5, 0.5, 0.5].into();
        assert!(a.point_hit(&p));
    }
}
