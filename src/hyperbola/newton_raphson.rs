use crate::hyperbola::util::{arsinh, sinh_cosh};

const MAX_ITERATIONS: usize = 20;

/// See https://arxiv.org/pdf/1503.01641v1
fn starter(eccentricity: f64, mean_anomaly: f64) -> f64 {
    let g = 1.0 / eccentricity;
    let l = mean_anomaly / eccentricity;
    let sinh_eccentric_anomaly = if 4.0 - 1.9*g < l {
        l + 2.3*g
    } else if 2.74 - 1.56*g < l {
        l + 1.9*g
    } else if 2.01 - 1.33*g < l {
        l + 1.56*g
    } else if 1.60 - 1.16*g < l {
        l + 1.33*g
    } else if 1.32 - 1.02*g < l {
        l + 1.16*g
    } else if 1.12 - 0.91*g < l {
        l + 1.02*g
    } else if 1.0 - 0.833333*g < l {
        l + 0.91*g
    } else {
        let a = f64::cbrt(3.0*l / g + f64::sqrt(9.0*l.powi(2) / g.powi(2) + 8.0*(1.0-g).powi(3) / g.powi(3)));
        a - 2.0*(1.0 - g) / (g * a)
    };
    arsinh(sinh_eccentric_anomaly)
}

/// Simple Newton-Raphson solver. Reasonably fast for low accuracies, terribly slow for high
/// accuracies.
///
/// Initial estimate comes from https://arxiv.org/pdf/1503.01641v1
pub fn solve_hke_with_newton_raphson(eccentricity: f64, mean_anomaly: f64, tolerance: f64) -> f64 {
    let mh = mean_anomaly.abs();
    let mut eccentric_anomaly = starter(eccentricity, mean_anomaly);

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

