fn sinh_cosh(a: f64) -> (f64, f64) {
    let ex = f64::exp(a);
    let enx = f64::exp(-a);
    let sinh = (ex-enx) / 2.0;
    let cosh = (ex+enx) / 2.0;
    (sinh, cosh)
}

pub fn hyperbolic_kepler_equation(eccentricity: f64, eccentric_anomaly: f64) -> f64 {
    eccentricity * f64::sinh(eccentric_anomaly) - eccentric_anomaly
}

pub fn solve_hyperbolic_kepler_equation(eccentricity: f64, mean_anomaly: f64) -> f64 {
    let mh = mean_anomaly.abs();
    let x = mh / eccentricity;
    let arsinh_x = f64::ln(x + f64::sqrt(x.powi(2) + 1.0));
    let mut eccentric_anomaly = arsinh_x;

    // Initial Halley iteration
    let (sinh_eccentric_anomaly, cosh_eccentric_anomaly) = sinh_cosh(eccentric_anomaly);
    let f = eccentricity * sinh_eccentric_anomaly - eccentric_anomaly - mh;
    let f_prime = eccentricity * cosh_eccentric_anomaly - 1.0;
    let f_prime_prime = f_prime + 1.0;
    eccentric_anomaly = eccentric_anomaly - (2.0 * f / f_prime) / (2.0 - f * f_prime_prime / f_prime.powi(2));

    // Newton-Raphson
    for i in 0..50 {
        let (sf0, cf0) = sinh_cosh(eccentric_anomaly);
        let error = eccentricity * sf0 - eccentric_anomaly - mh;
        if error.abs() < 1.0e-10 {
            break;
        }
        let error_prime = eccentricity * cf0 - 1.0;
        eccentric_anomaly = eccentric_anomaly - error / error_prime;
    }

    eccentric_anomaly
}

#[cfg(test)]
mod test {
    use crate::hyperbola::{hyperbolic_kepler_equation, solve_hyperbolic_kepler_equation};

    #[test]
    fn test_hyperbola() {
        let eccentricities: Vec<f64> = (1..999)
            .map(|x| 1.0 + f64::powi(x as f64, 2) / 1000.0)
            .collect();
        let mean_anomalies: Vec<f64> = (0..10000)
            .map(|x| f64::powi(x as f64, 2) / 10000.0)
            .collect();

        for e in &eccentricities {
            for m in &mean_anomalies {
                let eccentric_anomaly = solve_hyperbolic_kepler_equation(*e, *m);
                let predicted_m = hyperbolic_kepler_equation(*e, eccentric_anomaly);
                let difference = (m - predicted_m).abs();
                if difference > 1.0e-4 {
                    dbg!(m, predicted_m, difference, e);
                    panic!()
                }
            }
        }
    }
}
