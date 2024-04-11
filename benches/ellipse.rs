use std::time::Duration;

use criterion::{criterion_group, criterion_main, Bencher, Criterion};
use rust_kepler_solver::ellipse::EllipseSolver;

pub fn bench(c: &mut Criterion) {
    let eccentricities = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99];
    let mean_anomalies = [0.01, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 6.283];

    let mut group = c.benchmark_group("ellipse");

    group.warm_up_time(Duration::from_millis(1000));
    group.measurement_time(Duration::from_millis(2000));

    for e in eccentricities {
        let solver = EllipseSolver::new(e);
        group.throughput(criterion::Throughput::Elements(mean_anomalies.len() as u64));
        group.bench_function(format!("{}", e).as_str(), |b: &mut Bencher| {
            b.iter(|| {
                for m in &mean_anomalies {
                    solver.solve(*m);
                };
            });
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);