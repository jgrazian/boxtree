# Rust-BVH

A basic **B**ounding **V**olume **H**ierarchy implementation using [glam](https://github.com/bitshifter/glam-rs).

Rust-BVH aims to be easy to use for any arbitrary 2d or 3d objects.

```toml
[dependencies]
rust-bvh = { version = "0.1.0" }
```

## Basic Example

First, implement `Bounded`.

```rust
use rust-bvh::{Bounded, Bounds2, Bvh2, Ray2, Vec2};

struct Circle {
    center: Vec2,
    radius: f32,
}

impl Bounded<Bounds2> for Circle {
    type Item = Self;
    
    fn bounds(&self) -> Bounds2 {
        Bounds2 {
            min: self.center - self.radius,
            max: self.center + self.radius,
        }
    }
}
```

Next build a BVH.

```rust
let circles = vec![
    Circle {
        center: Vec2::new(0.0, 0.0),
        radius: 1.0,
    },
    Circle {
        center: Vec2::new(2.0, 0.0),
        radius: 1.0,
    },
    Circle {
        center: Vec2::new(4.0, 0.0),
        radius: 1.0,
    },
];

let circle_bvh = Bvh2::build(circles);
```

Finally, query the BVH.

> Note: query_ functions return all objects contained in bounding boxes that pass the query. For stronger filtering of only objects that are actually hit extra trait impls are required (see below).

```rust
// Using a ray
let query_ray = Ray2::new(Vec2::new(0.0, -2.0), Vec2::new(0.0, 1.0));
let hits_iter = circle_bvh.query_ray(&query_ray, 0.0, f32::MAX);
let hit_circle = hits_iter.next().unwrap().1;
assert_eq!(hit_circle.center, Vec2::new(0.0, 0.0));

// Using a bounding box
let query_bounds = Bounds2::new(Vec2::new(0.0, 0.0), Vec2::new(1.0, 1.0));
let hits_iter = circle_bvh.query_bounds(&query_bounds);
let hit_circle = hits_iter.next().unwrap().1;
assert_eq!(hit_circle.center, Vec2::new(0.0, 0.0));
```

## Advanced Usage

Rust-BVH has specialized traits for fast traversal.

```rust
impl RayHittable<Bounds2> for Circle {
    fn ray_hit(&self, ray: &Bounds2::Ray, t_min: f32, t_max: f32) -> Option<(f32, &Self::Item)> {
        let om = self.center - ray.origin;
        let a = ray.direction.length_squared();
        let b = 2.0*ray.dot(om);
        let c = om.length_squared() - (self.radius * self.radius);
        
        let q = b*b - 4.0*a*c;
        let g = 1.0 / (2.0*a);
        let q = g * q.sqrt();
        let b = -b*g;
        
        let t0 = D*B + Q) + ray.origin;
        let t1 = D*(B - Q) + ray.origin;
    
    	t0.min(t1)
    }
} 

let query_ray = Ray2::new(Vec2::new(0.0, -2.0), Vec2::new(0.0, 1.0));
let (time, hit_circle) = circle_bvh.ray_hit(query_ray, 0.0, f32::MAX).unwrap(); // No iter this time we get an object directly
assert_eq!(hit_circle.center, Vec2::new(0.0, 0.0));
```
