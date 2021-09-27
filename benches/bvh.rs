use criterion::black_box;
use criterion::*;
use rust_bvh::*;

fn setup_bvh2() -> Bvh2<Bounds2> {
    let objects: Vec<Bounds2> = (0..1000)
        .map(|i| ([i as f32; 2], [(i + 1) as f32; 2]).into())
        .collect();

    Bvh2::build(objects)
}

fn setup_bvh3() -> Bvh3<Bounds3> {
    let objects: Vec<Bounds3> = (0..1000)
        .map(|i| ([i as f32; 3], [(i + 1) as f32; 3]).into())
        .collect();

    Bvh3::build(objects)
}

fn setup_bvh3a() -> Bvh3A<Bounds3A> {
    let objects: Vec<Bounds3A> = (0..1000)
        .map(|i| ([i as f32; 3], [(i + 1) as f32; 3]).into())
        .collect();

    Bvh3A::build(objects)
}

fn bvh_build(c: &mut Criterion) {
    let objects_2: Vec<Bounds2> = (0..1000)
        .map(|i| ([i as f32; 2], [(i + 1) as f32; 2]).into())
        .collect();
    let objects_3: Vec<Bounds3> = (0..1000)
        .map(|i| ([i as f32; 3], [(i + 1) as f32; 3]).into())
        .collect();
    let objects_3a: Vec<Bounds3A> = (0..1000)
        .map(|i| ([i as f32; 3], [(i + 1) as f32; 3]).into())
        .collect();

    c.bench_function("Bvh2 build", move |b| {
        // This will avoid timing the to_vec call.
        b.iter_batched(
            || objects_2.clone(),
            |data| Bvh2::build(data),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("Bvh3 build", move |b| {
        // This will avoid timing the to_vec call.
        b.iter_batched(
            || objects_3.clone(),
            |data| Bvh3::build(data),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("Bvh3A build", move |b| {
        // This will avoid timing the to_vec call.
        b.iter_batched(
            || objects_3a.clone(),
            |data| Bvh3A::build(data),
            BatchSize::SmallInput,
        )
    });
}

fn bvh_iter(c: &mut Criterion) {
    let bvh2 = setup_bvh2();
    let bvh3 = setup_bvh3();
    let bvh3a = setup_bvh3a();

    c.bench_function("Bvh2 iter", move |b| {
        b.iter(|| bvh2.iter_objects().map(|obj| black_box(obj)));
    });
    c.bench_function("Bvh3 iter", move |b| {
        b.iter(|| bvh3.iter_objects().map(|obj| black_box(obj)));
    });
    c.bench_function("Bvh3A iter", move |b| {
        b.iter(|| bvh3a.iter_objects().map(|obj| black_box(obj)));
    });
}

fn bvh_ray(c: &mut Criterion) {
    let bvh2 = setup_bvh2();
    let bvh3 = setup_bvh3();
    let bvh3a = setup_bvh3a();

    let ray2 = ([0.0; 2], [1.0, 0.9]).into();
    let ray3 = ([0.0; 3], [1.0, 1.0, 0.9]).into();
    let ray3a = ([0.0; 3], [1.0, 1.0, 0.9]).into();

    c.bench_function("Bvh2 query_ray_exact", move |b| {
        b.iter(|| {
            bvh2.query_ray_exact(&ray2, 0.0, f32::MAX)
                .map(|obj| black_box(obj))
        });
    });
    c.bench_function("Bvh3 query_ray_exact", move |b| {
        b.iter(|| {
            bvh3.query_ray_exact(&ray3, 0.0, f32::MAX)
                .map(|obj| black_box(obj))
        });
    });
    c.bench_function("Bvh3A query_ray_exact", move |b| {
        b.iter(|| {
            bvh3a
                .query_ray_exact(&ray3a, 0.0, f32::MAX)
                .map(|obj| black_box(obj))
        });
    });
}

fn bvh_bounds(c: &mut Criterion) {
    let bvh2 = setup_bvh2();
    let bvh3 = setup_bvh3();
    let bvh3a = setup_bvh3a();

    let b2 = ([0.0; 2], [100.0; 2]).into();
    let b3 = ([0.0; 3], [100.0; 3]).into();
    let b3a = ([0.0; 3], [100.0; 3]).into();

    c.bench_function("Bvh2 query_bounds_exact", move |b| {
        b.iter(|| bvh2.query_bounds_exact(&b2).map(|obj| black_box(obj)));
    });
    c.bench_function("Bvh3 query_bounds_exact", move |b| {
        b.iter(|| bvh3.query_bounds_exact(&b3).map(|obj| black_box(obj)));
    });
    c.bench_function("Bvh3A query_bounds_exact", move |b| {
        b.iter(|| bvh3a.query_bounds_exact(&b3a).map(|obj| black_box(obj)));
    });
}

fn bvh_point(c: &mut Criterion) {
    let bvh2 = setup_bvh2();
    let bvh3 = setup_bvh3();
    let bvh3a = setup_bvh3a();

    let p2 = [100.5; 2].into();
    let p3 = [100.5; 3].into();
    let p3a = [100.5; 3].into();

    c.bench_function("Bvh2 query_point_exact", move |b| {
        b.iter(|| bvh2.query_point_exact(&p2).map(|obj| black_box(obj)));
    });
    c.bench_function("Bvh3 query_point_exact", move |b| {
        b.iter(|| bvh3.query_point_exact(&p3).map(|obj| black_box(obj)));
    });
    c.bench_function("Bvh3A query_point_exact", move |b| {
        b.iter(|| bvh3a.query_point_exact(&p3a).map(|obj| black_box(obj)));
    });
}

criterion_group!(benches, bvh_build, bvh_iter, bvh_ray, bvh_bounds, bvh_point);
criterion_main!(benches);
