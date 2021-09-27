use crate::traits::*;
use crate::{Ray2, Ray3, Ray3A};

use glam::{Vec2, Vec3, Vec3A};

pub type Bounds2 = Bounds<Vec2>;
pub type Bounds3 = Bounds<Vec3>;
pub type Bounds3A = Bounds<Vec3A>;

/// A bounding box in D dimensional space.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Bounds<V> {
    pub min: V,
    pub max: V,
}

impl<V> Bounds<V> {
    pub fn new(min: V, max: V) -> Bounds<V> {
        Bounds { min, max }
    }
}

impl<V> From<(V, V)> for Bounds<V> {
    fn from(tuple: (V, V)) -> Self {
        Self::new(tuple.0, tuple.1)
    }
}

// ================================================
// Bounds2
// ================================================
impl BoundingBox for Bounds2 {
    const DIM: usize = 2;
    type Vector = Vec2;
    type Ray = Ray2;

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
        shape.x * shape.y
    }

    fn union(&self, other: &Self) -> Self {
        Self {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }
}

impl Bounded<Bounds2> for Bounds2 {
    type Item = Bounds2;

    fn bounds(&self) -> Bounds2 {
        *self
    }
}

impl From<([f32; 2], [f32; 2])> for Bounds2 {
    fn from(tuple: ([f32; 2], [f32; 2])) -> Self {
        Self::new(tuple.0.into(), tuple.1.into())
    }
}

impl RayHittable<Bounds2> for Bounds2 {
    fn ray_hit(
        &self,
        ray: &<Self as BoundingBox>::Ray,
        t_min: f32,
        t_max: f32,
    ) -> Option<(f32, &Self::Item)> {
        let inv_d = ray.direction.recip();
        let t0__ = (self.min - ray.origin) * inv_d;
        let t1__ = (self.max - ray.origin) * inv_d;

        let mask = inv_d.cmplt(<Self as BoundingBox>::Vector::ZERO);
        let t0_ = <Self as BoundingBox>::Vector::select(mask, t1__, t0__);
        let t1_ = <Self as BoundingBox>::Vector::select(mask, t0__, t1__);

        let t0 = t0_.max_element();
        let t1 = t1_.min_element();

        dbg!(t0, t1);

        if t1 <= t0 {
            return None;
        }

        let c1 = t_min <= t0;
        let c2 = t_min <= t1;
        let c3 = t_max >= t0;
        let c4 = t_max >= t1;

        match (c1, c2, c3, c4) {
            (true, _, false, _) => None,
            (true, _, true, _) => Some((t0, &self)),
            (false, true, _, false) => None,
            (false, true, _, true) => Some((t1, &self)),
            (_, false, _, _) => None,
        }
    }
}

impl BoundsHittable<Bounds2> for Bounds2 {
    fn bounds_hit(&self, bounds: &Self) -> bool {
        (self.min.cmple(bounds.max) & self.max.cmpge(bounds.min)).all()
    }
}

impl PointHittable<Bounds2> for Bounds2 {
    fn point_hit(&self, point: &<Self as BoundingBox>::Vector) -> bool {
        (self.min.cmple(*point) & self.max.cmpge(*point)).all()
    }
}

// ================================================
// Bounds3
// ================================================
impl BoundingBox for Bounds3 {
    const DIM: usize = 3;
    type Vector = Vec3;
    type Ray = Ray3;

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
        2.0 * (shape.x * shape.y + shape.x * shape.z + shape.y * shape.z)
    }

    fn union(&self, other: &Self) -> Self {
        Self {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }
}

impl Bounded<Bounds3> for Bounds3 {
    type Item = Bounds3;

    fn bounds(&self) -> Bounds3 {
        *self
    }
}

impl From<([f32; 3], [f32; 3])> for Bounds3 {
    fn from(tuple: ([f32; 3], [f32; 3])) -> Self {
        Self::new(tuple.0.into(), tuple.1.into())
    }
}

impl RayHittable<Bounds3> for Bounds3 {
    fn ray_hit(
        &self,
        ray: &<Self as BoundingBox>::Ray,
        t_min: f32,
        t_max: f32,
    ) -> Option<(f32, &Self::Item)> {
        let inv_d = ray.direction.recip();
        let t0__ = (self.min - ray.origin) * inv_d;
        let t1__ = (self.max - ray.origin) * inv_d;

        let mask = inv_d.cmplt(<Self as BoundingBox>::Vector::ZERO);
        let t0_ = <Self as BoundingBox>::Vector::select(mask, t1__, t0__);
        let t1_ = <Self as BoundingBox>::Vector::select(mask, t0__, t1__);

        let t0 = t0_.max_element();
        let t1 = t1_.min_element();

        if t1 <= t0 {
            return None;
        }

        let c1 = t_min <= t0;
        let c2 = t_min <= t1;
        let c3 = t_max >= t0;
        let c4 = t_max >= t1;

        match (c1, c2, c3, c4) {
            (true, _, false, _) => None,
            (true, _, true, _) => Some((t0, &self)),
            (false, true, _, false) => None,
            (false, true, _, true) => Some((t1, &self)),
            (_, false, _, _) => None,
        }
    }
}

impl BoundsHittable<Bounds3> for Bounds3 {
    fn bounds_hit(&self, bounds: &Self) -> bool {
        (self.min.cmple(bounds.max) & self.max.cmpge(bounds.min)).all()
    }
}

impl PointHittable<Bounds3> for Bounds3 {
    fn point_hit(&self, point: &<Self as BoundingBox>::Vector) -> bool {
        (self.min.cmple(*point) & self.max.cmpge(*point)).all()
    }
}

// ================================================
// Bounds3A
// ================================================
impl BoundingBox for Bounds3A {
    const DIM: usize = 3;
    type Vector = Vec3A;
    type Ray = Ray3A;

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
        2.0 * (shape.x * shape.y + shape.x * shape.z + shape.y * shape.z)
    }

    fn union(&self, other: &Self) -> Self {
        Self {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }
}

impl Bounded<Bounds3A> for Bounds3A {
    type Item = Bounds3A;

    fn bounds(&self) -> Bounds3A {
        *self
    }
}

impl From<([f32; 3], [f32; 3])> for Bounds3A {
    fn from(tuple: ([f32; 3], [f32; 3])) -> Self {
        Self::new(tuple.0.into(), tuple.1.into())
    }
}

impl RayHittable<Bounds3A> for Bounds3A {
    fn ray_hit(
        &self,
        ray: &<Self as BoundingBox>::Ray,
        t_min: f32,
        t_max: f32,
    ) -> Option<(f32, &Self::Item)> {
        let inv_d = ray.direction.recip();
        let t0__ = (self.min - ray.origin) * inv_d;
        let t1__ = (self.max - ray.origin) * inv_d;

        let mask = inv_d.cmplt(<Self as BoundingBox>::Vector::ZERO);
        let t0_ = <Self as BoundingBox>::Vector::select(mask, t1__, t0__);
        let t1_ = <Self as BoundingBox>::Vector::select(mask, t0__, t1__);

        let t0 = t0_.max_element();
        let t1 = t1_.min_element();

        if t1 <= t0 {
            return None;
        }

        let c1 = t_min <= t0;
        let c2 = t_min <= t1;
        let c3 = t_max >= t0;
        let c4 = t_max >= t1;

        match (c1, c2, c3, c4) {
            (true, _, false, _) => None,
            (true, _, true, _) => Some((t0, &self)),
            (false, true, _, false) => None,
            (false, true, _, true) => Some((t1, &self)),
            (_, false, _, _) => None,
        }
    }
}

impl BoundsHittable<Bounds3A> for Bounds3A {
    fn bounds_hit(&self, bounds: &Self) -> bool {
        (self.min.cmple(bounds.max) & self.max.cmpge(bounds.min)).all()
    }
}

impl PointHittable<Bounds3A> for Bounds3A {
    fn point_hit(&self, point: &<Self as BoundingBox>::Vector) -> bool {
        (self.min.cmple(*point) & self.max.cmpge(*point)).all()
    }
}

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
        let r1 = ([0.5, 2.0], [0.0, -1.0]).into();
        let r2 = ([0.5, 0.5], [0.0, 1.0]).into();
        let r3 = ([0.5, 2.0], [1.0, 0.0]).into();
        let r4 = ([0.5, 2.0], [0.0, 1.0]).into();
        assert_eq!(a.ray_hit(&r1, 0.0, f32::MAX), Some((1.0, &a))); // hit
        assert_eq!(a.ray_hit(&r2, 0.0, f32::MAX), Some((0.5, &a))); // hit inside
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((-0.5, &a)));
        assert_eq!(a.ray_hit(&r3, 0.0, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, 0.0, f32::MAX), None);

        let r = ([0.0, -0.5], [0.0, 1.0]).into();
        assert_eq!(a.ray_hit(&r, 0.0, 0.4), None);
        assert_eq!(a.ray_hit(&r, 1.6, f32::MAX), None);
        assert_eq!(a.ray_hit(&r, 1.4, f32::MAX), Some((1.5, &a)));
    }

    #[test]
    fn ray_hit_3d() {
        let a: Bounds3 = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).into();
        let r1 = ([0.5, 2.0, 0.5], [0.0, -1.0, 0.0]).into(); // intersects
        let r2 = ([0.5, 0.5, 0.5], [0.0, 1.0, 0.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0, 0.5], [1.0, 0.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0, 0.5], [0.0, 1.0, 0.0]).into(); // miss
        assert_eq!(a.ray_hit(&r1, 0.0, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.ray_hit(&r2, 0.0, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((-0.5, &a)));
        assert_eq!(a.ray_hit(&r3, 0.0, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, 0.0, f32::MAX), None);

        let r = ([0.5, 0.5, -0.5], [0.0, 0.0, 1.0]).into();
        assert_eq!(a.ray_hit(&r, 0.0, 0.4), None);
        assert_eq!(a.ray_hit(&r, 1.6, f32::MAX), None);
        assert_eq!(a.ray_hit(&r, 1.4, f32::MAX), Some((1.5, &a)));
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
