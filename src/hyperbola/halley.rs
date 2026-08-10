use crate::hyperbola::util::sinh_cosh;

const MAX_ITERATIONS: usize = 200;

/// Simple Halley solver. Reasonably fast for low accuracies, terribly slow for high
/// accuracies.
///
/// Initial estimate comes from https://arxiv.org/pdf/1503.01641v1
pub fn solve_hke_with_halley(eccentricity: f64, mean_anomaly: f64, tolerance: f64) -> f64 {
    let mh = mean_anomaly.abs();
    let x = mh / eccentricity;
    let arsinh_x = f64::ln(x + f64::sqrt(x.powi(2) + 1.0));
    let mut eccentric_anomaly = arsinh_x;

    // Initial Halley iteration
    //let (sinh_eccentric_anomaly, cosh_eccentric_anomaly) = sinh_cosh(eccentric_anomaly);
    //let f = eccentricity * sinh_eccentric_anomaly - eccentric_anomaly - mh;
    //let f_prime = eccentricity * cosh_eccentric_anomaly - 1.0;
    //let f_prime_prime = f_prime + 1.0;
    //eccentric_anomaly = eccentric_anomaly - (2.0 * f / f_prime) / (2.0 - f * f_prime_prime / f_prime.powi(2));

    for _ in 0..MAX_ITERATIONS {
        let (sf0, cf0) = sinh_cosh(eccentric_anomaly);
        let error = eccentricity * sf0 - eccentric_anomaly - mh;
        if error.abs() < tolerance {
            break;
        }
        let error_prime = eccentricity * cf0 - 1.0;
        eccentric_anomaly = eccentric_anomaly - error / error_prime;
    }

    eccentric_anomaly
}

