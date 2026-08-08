use lazy_static::lazy_static;
use serde::{Deserialize, Serialize};

const CUBIC_DELTA_THRESHOLD: f32 = 1.0e-6;

lazy_static! {
    // From eq. 4 in the B. Wu et all paper
    // Max value in each interval
    static ref PADE_ECCENTRIC_ANOMALY_THRESHOLDS: [f32; 15] = [
        40.0 / 8.0,
        38.0 / 8.0,
        34.0 / 8.0,
        30.0 / 8.0,
        26.0 / 8.0,
        22.0 / 8.0,
        18.0 / 8.0,
        15.0 / 8.0,
        13.0 / 8.0,
        11.0 / 8.0,
        9.0  / 8.0,
        7.0  / 8.0,
        5.0  / 8.0,
        3.0  / 8.0,
        29.0 / 200.0
    ];

    static ref PADE_ORDERS: [f32; 15] = [
        10.0 / 2.0,
        9.0 / 2.0,
        8.0 / 2.0,
        7.0 / 2.0,
        6.0 / 2.0,
        5.0 / 2.0,
        8.0 / 4.0,
        7.0 / 4.0,
        6.0 / 4.0,
        5.0 / 4.0,
        4.0 / 4.0,
        3.0 / 4.0,
        2.0 / 4.0,
        1.0 / 4.0,
        0.0 / 4.0
    ];
}

fn sinh_cosh(a: f32) -> (f32, f32) {
    let a = a as f64;
    let ex = f64::exp(a);
    let enx = f64::exp(-a);
    let sinh = (ex-enx) as f32 / 2.0;
    let cosh = (ex+enx) as f32 / 2.0;
    (sinh, cosh)
}

fn fast_sinh_cosh(a: f32) -> (f32, f32) {
    let ex = fast_math::exp(a);
    let enx = fast_math::exp(-a);
    let sinh = (ex-enx) / 2.0;
    let cosh = (ex+enx) / 2.0;
    (sinh, cosh)
}

fn hyperbolic_kepler_equation(eccentricity: f32, eccentric_anomaly: f32) -> f32 {
    eccentricity * f32::sinh(eccentric_anomaly) - eccentric_anomaly
}

fn pade_approximation(ec: f32, mh: f32, a: f32) -> [f32; 4] {
    let (sa, ca) = fast_sinh_cosh(a);
    let d1 = ca.powi(2) + 3.0;
    let d2 = sa.powi(2) + 4.0;
    let p1 = ca * (3.0 * ca.powi(2) + 17.0) / (5.0 * d1);
    let p2 = sa * (3.0 * sa.powi(2) + 28.0) / (20.0 * d2);
    let p3 = ca * (ca.powi(2) + 27.0) / (60.0 * d1);
    let q1 = -2.0 * ca * sa / (5.0 * d1);
    let q2 = (sa.powi(2) - 4.0) / (20.0 * d2);
    let coefficient_f3 = ec * p3 - q2;
    let coefficient_f2 = ec * p2 - (mh + a) * q2 - q1;
    let coefficient_f1 = ec * p1 - (mh + a) * q1 - 1.0;
    let coefficient_f0 = ec * sa - mh - a;
    [coefficient_f3, coefficient_f2, coefficient_f1, coefficient_f0]
}

/// f(x) = cf[0]*x^3 + cf[1]*x^2 + cf[2]*x + cf[3]
fn solve_cubic(coefficients: [f32; 4], mh: f32, ec: f32) -> f32 {
    let mut x = f32::min(mh, ec) / (ec - 0.538);
    loop {
        // halley's method
        let f = ((coefficients[0]*x + coefficients[1])*x + coefficients[2])*x + coefficients[3];
        let f_prime = (3.0*coefficients[0]*x + 2.0*coefficients[1])*x + coefficients[2];
        let f_prime_prime = 6.0*coefficients[0]*x + 2.0*coefficients[1];
        let delta = -2.0*f*f_prime / (2.0*f_prime.powi(2) - f*f_prime_prime);
        if delta < CUBIC_DELTA_THRESHOLD && delta > -CUBIC_DELTA_THRESHOLD {
            break;
        }
        x += delta;
    }
    x
}

/// ## Example
/// ```rs
/// use rust_kepler_solver::hyperbola::HyperbolaSolver;
///
/// fn example_hyperbola() {
///     let eccentricity = 1.0;
///     let solver = HyperbolaSolver::new(eccentricity);
///     println!("{}", solver.solve(1.2));
///     println!("{}", solver.solve(100.0));
/// }
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HyperbolaSolver {
    eccentricity: f32,
    pade_mean_anomaly_thresholds: Vec<f32>,
}

impl HyperbolaSolver {
    pub fn new(eccentricity: f32) -> Self {
        assert!(eccentricity >= 1.0);
        // Note: the mean anomaly thresholds are not actually monotonically increasing and for
        // lower values of e seem to curve back upwards towards the end. I am not sure why this is
        // and it does not seem to be addressed in the paper, but regardless the approximation
        // seems to be very good so we deal with this by just terminating if the values start to
        // curve back up.
        let mut pade_mean_anomaly_thresholds = vec![];
        for eccentric_anomaly in PADE_ECCENTRIC_ANOMALY_THRESHOLDS.iter() {
            let mean_anomaly = hyperbolic_kepler_equation(eccentricity, *eccentric_anomaly);
            if let Some(previous_mean_anomaly) = pade_mean_anomaly_thresholds.last() && mean_anomaly > *previous_mean_anomaly {
                break;
            }
            pade_mean_anomaly_thresholds.push(mean_anomaly);
        }

        Self { eccentricity, pade_mean_anomaly_thresholds: pade_mean_anomaly_thresholds.into() }
    }

    fn get_initial_approximation_for_small_mh(&self, ec: f32, mh: f32) -> f32 {
        // For mh < hke(5) we use a 'piecewise pade approximation' to get the starting estimate
        let i = self.pade_mean_anomaly_thresholds.partition_point(|&x| x >= mh)-1;
        let a = PADE_ORDERS[i];
        let coefficients = pade_approximation(ec, mh, a);

        solve_cubic(coefficients, mh, ec) + a
    }

    fn get_initial_approximation_for_large_mh(&self, ec: f32, mh: f32) -> f32 {
        // For mh >= hke(5), we can use this... thing that I copied from the above paper
        // I have no idea how it works, but it works very very very well
        let fa = f32::ln(2.0 * mh / ec);
        let ca = 0.5 * (2.0 * mh / ec + ec / (2.0 * mh));
        let sa = 0.5 * (2.0 * mh / ec - ec / (2.0 * mh));
        let top = 6.0 * (ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0) 
        + 3.0 * (ec * sa / (ec * ca - 1.0)) * ((ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0)).powi(2);
        let bottom = 6.0 + 6.0 * (ec * sa / (ec * ca - 1.0)) * ((ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0))
        + (ec * ca / (ec * ca - 1.0)) * ((ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0)).powi(2);
        let delta = top / bottom;
        fa + delta
    }

    fn improve_initial_approximation(&self, f0: f32, mh: f32, ec: f32, mean_anomaly: f32) -> f32 {
        // Halley method
        let (sf0, cf0) = sinh_cosh(f0);
        let f = ec * sf0 - f0 - mh;
        let f_prime = ec * cf0 - 1.0;
        let f_prime_prime = f_prime + 1.0;
        let f1 = f0 - (2.0 * f / f_prime) / (2.0 - f * f_prime_prime / f_prime.powi(2));
        f1 * mean_anomaly.signum()
    }

    /// Works with all values of mean anomaly 0 to infinity
    pub fn solve(&self, mean_anomaly: f32) -> f32 {
        // Solver assumes mean anomaly > 0
        // The equation is symmetric, so for mean anomaly < 0, we just flip the sign o the output
        let ec = self.eccentricity;
        let mh = mean_anomaly.abs();

        let f0 = if mh <= self.pade_mean_anomaly_thresholds[0] {
            self.get_initial_approximation_for_small_mh(ec, mh)
        } else {
            self.get_initial_approximation_for_large_mh(ec, mh)
        };
        self.improve_initial_approximation(f0, mh, ec, mean_anomaly)
    }
}

#[cfg(test)]
mod test {
    use crate::bisection::bisection;

    use super::HyperbolaSolver;

    fn solve_with_bisection(e: f32, m: f32) -> f32 {
        // We don't care about speed here, so just use as wide a range as possible
        let f = |eccentric_anomaly: f32| e * f32::sinh(eccentric_anomaly) - eccentric_anomaly - m;
        bisection(&f, -100000.0, 100000.0)
    }

    #[test]
    fn test_hyperbola() {
        let eccentricities: Vec<f32> = (1..999)
            .map(|x| 1.0 + f32::powi(x as f32, 2) / 1000.0)
            .collect();
        let mean_anomalies: Vec<f32> = (0..10000)
            .map(|x| f32::powi(x as f32, 2) / 10000.0)
            .collect();
            
        for e in &eccentricities {
            let solver = HyperbolaSolver::new(*e);
            for m in &mean_anomalies {
                let expected = solve_with_bisection(*e, *m);
                let actual = solver.solve(*m);
                let difference = if actual.abs() == 0.0 { expected - actual } else { (expected - actual) / actual }.abs();
                if difference > 1.0e-4 {
                    dbg!(expected, actual, difference, e, m);
                    panic!()
                }
            }
        }
    }
}
