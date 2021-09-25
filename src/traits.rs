pub trait Bounded<V>: Clone {
    type Item: Bounded<V>;
    type Bounds: Bounded<V>;

    fn bounds(&self) -> Self::Bounds;
}

pub trait RayHittable<V, R>: Bounded<V> {
    fn ray_hit(&self, ray: &R, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)>;
}

pub trait BoundsHittable<V>: Bounded<V> {
    fn bounds_hit(&self, bounds: &Self::Bounds) -> bool;
}

pub trait PointHittable<V>: Bounded<V> {
    fn point_hit(&self, point: &V) -> bool;
}
