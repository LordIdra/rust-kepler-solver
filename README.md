## Overview
This repository contains Rust solvers for the Elliptic Kepler Equation (EKE) and Hyperbolic Kepler Equation (HKE). Use `cargo test` in root to run tests, and `cargo bench` in root to run benchmarks. Note: This implementation uses f64's, greater performance can be achieved with f32's.

## Example
```rs
use std::f64::consts::PI;
use rust_kepler_solver::ellipse::EllipseSolver;

fn example_ellipse() {
    let eccentricity = 1.0;
    let solver = EllipseSolver::new(eccentricity);
    println!("{}", solver.solve(1.2));
    println!("{}", solver.solve(PI / 4.0));
}
```

```rs
use rust_kepler_solver::hyperbola::HyperbolaSolver;

fn example_hyperbola() {
    let eccentricity = 1.0;
    let solver = HyperbolaSolver::new(eccentricity);
    println!("{}", solver.solve(1.2));
    println!("{}", solver.solve(100.0));
}
```

## Method
### EKE
The EKE is solved by choosing an initial seed as described by Daniele Tommasini and David N. Olivieri (https://doi.org/10.1051/0004-6361/20214142), and then using Laguerre's method to iterate until the delta falls below a certain threshold. Laguerre's method is a reliable algorithm for solving the EKE according to Bruce A. Conway (https://doi.org/10.1007/BF01230852). There is almost certainly a more efficient method out there, but this implementation is still very fast.

### HKE
The HKE is solved with a slightly more complicated method as per Baisheng Wu et al (https://doi.org/10.1016/j.apm.2023.12.017). This method splits the interval of eccentric anomalies into two parts: one finite and one infinite part. An approximation is constructed for each region, the first using a piecewise Pade approximation, the second using 'an analytical initial approximate solution of the HKE.' We then compute thresholds for which interval a given mean anomaly should use, and get an initial approximation based off that. The approximations are so ridiculously accurate that only one step of Halley iteration is required to get a very precise result.

## Reliability
The crate includes tests for both the EKE and HKE solvers, which test ~ 6,000,000 and ~10,000,000 eccentricity and mean anomaly pairs. The values are linearly distributed for the EKE to cover the range of possible eccentricities and mean anomalies. For the HKE, both eccentricity and mean anomaly inputs up to infinity are technically valid, so we generate values using x^2/c to test a range of the smaller values (which is where the Pade approximation comes in) and larger values (where the analytical approximation comes in). Though it's not completely comprehensive, this should be enough to show that both solvers are very reliable.