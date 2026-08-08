use std::time::Duration;

use criterion::{Bencher, Criterion, criterion_group, criterion_main};
use rust_kepler_solver::hyperbola::HyperbolaSolver;

pub fn bench(c: &mut Criterion) {
    let eccentricities = [1.0001, 1.001, 1.01, 1.1, 11.0, 111.0, 1111.0];
    let mean_anomalies = [0.01, 0.5, 1.0, 2.0, 4.0, 8.0, 32.0, 128.0, 999.0, 9999999.0];

    let mut group = c.benchmark_group(format!("hyperbola"));

    group.warm_up_time(Duration::from_millis(20));
    group.measurement_time(Duration::from_millis(80));
    group.throughput(criterion::Throughput::Elements(mean_anomalies.len() as u64));

    for e in eccentricities {
        let solver = HyperbolaSolver::new(e);
        group.bench_function(format!("{}", e).as_str(), |b: &mut Bencher| {
            for m in mean_anomalies {
                b.iter(|| {
                    solver.solve(m);
                });
            }
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
