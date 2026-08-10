mod halley;
mod hkesdg;
mod newton_raphson;
mod util;

pub use hkesdg::solve_hke_with_hkesdg;
pub use newton_raphson::solve_hke_with_newton_raphson;

pub fn hke(eccentricity: f64, eccentric_anomaly: f64) -> f64 {
    eccentricity * f64::sinh(eccentric_anomaly) - eccentric_anomaly
}

#[cfg(test)]
mod test {
    use crate::hyperbola::{hke, solve_hke_with_hkesdg, solve_hke_with_newton_raphson};

    fn test_hke_solver(max_difference: f64, solver: impl Fn(f64, f64) -> f64) {
        let eccentricities: Vec<f64> = (1..999)
            .map(|x| 1.0 + f64::powi(x as f64, 2) / 1000.0)
            .collect();
        let mean_anomalies: Vec<f64> = (0..10000)
            .map(|x| f64::powi(x as f64, 2) / 10000.0)
            .collect();

        for e in &eccentricities {
            for m in &mean_anomalies {
                let eccentric_anomaly = solver(*e, *m);
                let predicted_m = hke(*e, eccentric_anomaly);
                let difference = (m - predicted_m).abs();
                if difference > max_difference {
                    dbg!(m, predicted_m, difference, e);
                    panic!()
                }
            }
        }
    }

    #[test]
    fn test_hke_hkesdg() {
        test_hke_solver(1.0e-5, |e: f64, m: f64| solve_hke_with_hkesdg(e, m, 1.0e-6));
        test_hke_solver(1.0e-10, |e: f64, m: f64| solve_hke_with_hkesdg(e, m, 1.0e-16));
    }

    #[test]
    fn test_hke_newton_raphson() {
        test_hke_solver(1.0e-5, |e: f64, m: f64| solve_hke_with_newton_raphson(e, m, 1.0e-6));
        test_hke_solver(1.0e-10, |e: f64, m: f64| solve_hke_with_newton_raphson(e, m, 1.0e-16));
    }
}
