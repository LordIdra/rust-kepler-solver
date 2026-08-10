use std::ffi::{c_double, c_int};

#[link(name = "hke-sdg")]
unsafe extern "C" {
    fn hkepler_eq(M1: c_double, e1: c_double, ph1: &mut c_double, resto1: &mut c_double, eps1: c_double) -> c_int;
}

/**
 * Invokes the [HKE-SDG solver](https://doi.org/10.1051/0004-6361/201833563) for the hyperbolic
 * kepler equation.
 *
 * Use the `tolerance` to tweak the accuracy of the result. This solver is comparatively slow 
 * for low accuracies, but can achieve extremely high accuracies in similar times.
 **/
pub fn solve_hke_with_hkesdg(
    eccentricity: f64, 
    mean_anomaly: f64, 
    tolerance: f64
) -> f64 {
    unsafe {
        let mut hyperbolic_anomaly = 0.0;
        let mut residual = 0.0;
        hkepler_eq(mean_anomaly, eccentricity, &mut hyperbolic_anomaly, &mut residual, tolerance);
        hyperbolic_anomaly
    }
}

