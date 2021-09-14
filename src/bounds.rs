use crate::traits::*;
use crate::{Ray, Vector};

/// A bounding box in D dimensional space.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Bounds<const D: usize> {
    pub min: Vector<D>,
    pub max: Vector<D>,
}

impl<const D: usize> Bounds<D> {
    pub fn new(min: Vector<D>, max: Vector<D>) -> Bounds<D> {
        Bounds { min, max }
    }

    pub fn shape(&self) -> Vector<D> {
        let mut shape = [0.0f32; D];
        (0..D).for_each(|i| shape[i] = self.max[i] - self.min[i]);
        shape
    }

    pub fn axis_length(&self, axis: usize) -> f32 {
        self.max[axis] - self.min[axis]
    }

    pub fn centroid(&self) -> Vector<D> {
        let mut centroid = [0.0; D];
        (0..D).for_each(|i| centroid[i] = (self.min[i] + self.max[i]) / 2.0);
        centroid
    }

    pub fn surface_area(&self) -> f32 {
        let mut total_area = 0.0;
        for i in 0..D {
            let length = self.axis_length(i);
            for j in (i + 1)..D {
                let length2 = self.axis_length(j);
                total_area += length * length2;
            }
        }
        (D - 1) as f32 * total_area
    }

    pub fn union(&self, other: &Bounds<D>) -> Bounds<D> {
        let mut min = [0.0; D];
        let mut max = [0.0; D];

        (0..D).for_each(|i| {
            min[i] = self.min[i].min(other.min[i]);
            max[i] = self.max[i].max(other.max[i])
        });

        Bounds { min: min, max: max }
    }
}

impl<const D: usize> Bounded<D> for Bounds<D> {
    type Bound = Bounds<D>;

    fn bounds(&self) -> Bounds<D> {
        *self
    }
}

impl<const D: usize> RayHittable<D> for Bounds<D> {
    fn ray_hit(&self, ray: &Ray<D>, t_min: f32, t_max: f32) -> Option<(f32, &Self::Bound)> {
        let mut loc_t_min = t_min;
        let mut loc_t_max = t_max;

        for a in 0..D {
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

impl<const D: usize> BoundsHittable<D> for Bounds<D> {
    fn bounds_hit(&self, bounds: &Bounds<D>) -> bool {
        (0..D).all(|i| self.min[i] <= bounds.max[i] && self.max[i] >= bounds.min[i])
    }
}

impl<const D: usize> PointHittable<D> for Bounds<D> {
    fn point_hit(&self, point: &Vector<D>) -> bool {
        (0..D).all(|i| self.min[i] <= point[i] && self.max[i] >= point[i])
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn axis_length_2d() {
        let b = Bounds::new([0.0, 0.0], [1.2, 1.3]);
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
    }

    #[test]
    fn axis_length_3d() {
        let b = Bounds::new([0.0, 0.0, 0.0], [1.2, 1.3, 1.4]);
        assert_eq!(b.axis_length(0), 1.2);
        assert_eq!(b.axis_length(1), 1.3);
        assert_eq!(b.axis_length(2), 1.4);
    }

    #[test]
    fn surface_area_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        assert_eq!(a.surface_area(), 1.0);
    }

    #[test]
    fn surface_area_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(a.surface_area(), 6.0);
    }

    #[test]
    fn union_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let b = Bounds::new([1.0, 1.0], [2.0, 2.0]);
        assert_eq!(a.union(&b), Bounds::new([0.0, 0.0], [2.0, 2.0]));
    }

    #[test]
    fn union_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Bounds::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        assert_eq!(a.union(&b), Bounds::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]));
    }

    #[test]
    fn ray_hit_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let r1 = ([0.5, 2.0], [0.0, -1.0]).into(); // intersects
        let r2 = ([0.5, 0.5], [0.0, 1.0]).into(); // intersects from inside
        let r3 = ([0.5, 2.0], [1.0, 0.0]).into(); // miss
        let r4 = ([0.5, 2.0], [0.0, 1.0]).into(); // miss
        assert_eq!(a.ray_hit(&r1, f32::MIN, f32::MAX), Some((1.0, &a)));
        assert_eq!(a.ray_hit(&r2, f32::MIN, f32::MAX), Some((0.5, &a)));
        assert_eq!(a.ray_hit(&r3, f32::MIN, f32::MAX), None);
        assert_eq!(a.ray_hit(&r4, f32::MIN, f32::MAX), None);

        let b = Bounds::new([3.0, -2.0], [4.0, -1.0]);
        let r = ([4.5, -2.5], [1.0, 1.0]).into();
        assert_eq!(b.ray_hit(&r, f32::MIN, f32::MAX), None);

        let b = Bounds::new([4.5, 0.0], [5.5, 1.0]);
        let r = ([4.5, -2.5], [1.0, 1.0]).into();
        assert_eq!(b.ray_hit(&r, f32::MIN, f32::MAX), None);
    }

    #[test]
    fn ray_hit_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
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
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let b = Bounds::new([2.0, 2.0], [3.0, 3.0]);
        assert!(!a.bounds_hit(&b));
        let c = Bounds::new([0.5, 0.5], [1.5, 1.5]);
        assert!(a.bounds_hit(&c));
    }

    #[test]
    fn bounds_hit_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Bounds::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(!a.bounds_hit(&b));
        let c = Bounds::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        assert!(a.bounds_hit(&c));
    }

    #[test]
    fn point_hit_2d() {
        let a = Bounds::new([0.0, 0.0], [1.0, 1.0]);
        let p = [2.0, 2.0];
        assert!(!a.point_hit(&p));
        let p = [0.5, 0.5];
        assert!(a.point_hit(&p));
    }

    #[test]
    fn point_hit_3d() {
        let a = Bounds::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let p = [2.0, 2.0, 2.0];
        assert!(!a.point_hit(&p));
        let p = [0.5, 0.5, 0.5];
        assert!(a.point_hit(&p));
    }
}
