use std::{hint::black_box, time::Instant};

use crate::{ellipse::EllipseSolver, hyperbola::HyperbolaSolver};

#[cfg(test)]
mod bisection;
mod ellipse;
mod hyperbola;

// https://stackoverflow.com/questions/44338311/rust-benchmark-optimized-out for why we use black_box
pub fn benchmark_ellipse() {
    let eccentricites: Vec<f64> = (1..999)
        .map(|x| x as f64 / 1000.0)
        .collect();
    let mean_anomalies: Vec<f64> = (0..6283)
        .map(|x| x as f64 / 1000.0)
        .collect();

    let now = Instant::now();
    for e in &eccentricites {
        let solver = EllipseSolver::new(*e);
        for m in &mean_anomalies {
            black_box(solver.solve(*m));
        }
    }
    let elapsed = now.elapsed();
    
    let iterations = (eccentricites.len() * mean_anomalies.len()) as u128;
    println!("Ellipse solver: {:.3?} total | {} iterations | {:.3}ns per iteration | {}/s", elapsed, 
        iterations, elapsed.as_nanos() / iterations, iterations as f64 / elapsed.as_secs_f64());
}

pub fn benchmark_hyperbola() {
    let eccentricites: Vec<f64> = (1..999)
        .map(|x| 1.0 + f64::powi(x as f64, 2) / 1000.0)
        .collect();
    let mean_anomalies: Vec<f64> = (0..10000)
        .map(|x| f64::powi(x as f64, 2) / 10000.0)
        .collect();

    let now = Instant::now();
    for e in &eccentricites {
        let solver = HyperbolaSolver::new(*e);
        for m in &mean_anomalies {
            black_box(solver.solve(*m));
        }
    }
    let elapsed = now.elapsed();

    let iterations = (eccentricites.len() * mean_anomalies.len()) as u128;
    println!("Hyperbola solver: {:.3?} total | {} iterations | {:.3}ns per iteration | {}/s", elapsed, 
        iterations, elapsed.as_nanos() / iterations, iterations as f64 / elapsed.as_secs_f64());
}

fn main() {
    benchmark_ellipse();
    benchmark_hyperbola();
}
