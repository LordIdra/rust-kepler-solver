use std::time::Duration;

use criterion::{BatchSize, Bencher, Criterion, Throughput, black_box, criterion_group, criterion_main};
use rand::seq::IndexedRandom;
use rust_kepler_solver::hyperbola::{solve_hke_with_hkesdg, solve_hke_with_newton_raphson};

pub fn bench(c: &mut Criterion) {
    let eccentricities = [1.0001, 1.01, 1.1, 11.0, 1111.0];
    let mean_anomalies = [0.01, 0.5, 1.0, 2.0, 4.0, 8.0, 32.0, 128.0, 999.0, 9999999.0];

    let mut group = c.benchmark_group(format!("hyperbola"));

    group.warm_up_time(Duration::from_millis(200));
    group.measurement_time(Duration::from_millis(800));
    group.throughput(Throughput::Elements(eccentricities.len() as u64));

    // HKE-SDG high accuracy
    for e in eccentricities {
        group.bench_function(format!("HKE-SDG (high accuracy) e={e}))"), |b: &mut Bencher| {
            b.iter_batched(
                || mean_anomalies.choose(&mut rand::rng()).unwrap(), 
                |&m| black_box(solve_hke_with_hkesdg(e, m, 1.0e-16)), 
                BatchSize::SmallInput
            );
        });
    }

    // HKE-SDG low accuracy
    for e in eccentricities {
        group.bench_function(format!("HKE-SDG (low accuracy) e={e}"), |b: &mut Bencher| {
            b.iter_batched(
                || mean_anomalies.choose(&mut rand::rng()).unwrap(), 
                |&m| black_box(solve_hke_with_hkesdg(e, m, 1.0e-5)), 
                BatchSize::SmallInput
            );
        });
    }

    // Newton-Raphson
    for e in eccentricities {
        group.bench_function(format!("Newton-Raphson (low accuracy) e={e}"), |b: &mut Bencher| {
            b.iter_batched(
                || mean_anomalies.choose(&mut rand::rng()).unwrap(), 
                |&m| black_box(solve_hke_with_newton_raphson(e, m, 1.0e-5)), 
                BatchSize::SmallInput
            );
        });
    }

    // Newton-Raphson
    for e in eccentricities {
        group.bench_function(format!("Newton-Raphson (medium accuracy) e={e}"), |b: &mut Bencher| {
            b.iter_batched(
                || mean_anomalies.choose(&mut rand::rng()).unwrap(), 
                |&m| black_box(solve_hke_with_newton_raphson(e, m, 1.0e-10)), 
                BatchSize::SmallInput
            );
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
