use lazy_static::lazy_static;

const CUBIC_DELTA_THRESHOLD: f64 = 1.0e-6;

lazy_static! {
    // From eq. 4 in the above paper
    // Max value in each interval
    static ref PADE_ECCENTRIC_ANOMALY_THRESHOLDS: [f64; 15] = [
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
        29.0 / 200.0];

    static ref PADE_ORDERS: [f64; 15] = [
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
        0.0 / 4.0];
}

fn hyperbolic_kepler_equation(eccentricity: f64, eccentric_anomaly: f64) -> f64 {
    eccentricity * f64::sinh(eccentric_anomaly) - eccentric_anomaly
}

fn pade_approximation(ec: f64, mh: f64, a: f64) -> [f64; 4] {
    let ex = f64::exp(a);
    let enx = f64::exp(-a);
    let sa = (ex-enx) / 2.0; // sinh(a)
    let ca = (ex+enx) / 2.0; // cosh(a)
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
    [coefficient_f3,coefficient_f2, coefficient_f1, coefficient_f0]
}

fn solve_cubic(coefficients: [f64; 4], mh: f64, ec: f64) -> f64 {
    let mut x = mh / (ec - 1.0); // starting value from series expansion of HKE
    loop {
        // halley's method
        let f = ((coefficients[0]*x + coefficients[1])*x + coefficients[2])*x + coefficients[3];
        let f_prime = (3.0*coefficients[0]*x + 2.0*coefficients[1])*x + coefficients[2];
        let f_prime_prime = 6.0*coefficients[0]*x + 2.0*coefficients[1];
        let delta = -2.0*f*f_prime / (2.0*f_prime.powi(2) - f*f_prime_prime);
        if delta.abs() < CUBIC_DELTA_THRESHOLD {
            break;
        }
        x += delta;
    }
    x
}

/// The method here comes from B. Wu et al, A new method for solving the hyperbolic Kepler equation
/// Basically, we divide up the range of eccentric anomalies into 0-5 and 5+
/// Then we compute what mean anomaly an eccentric anomaly of 5 comes out to, and cache that
/// If our mean anomaly is larger than this, we substitute into a gigantic expression
/// I don't really know how it works but it's VERY accurate for E > 5
/// For mean anomalies smaller than our threshold, we use a Pade approximation
#[derive(Debug)]
pub struct HyperbolaSolver {
    eccentricity: f64,
    pade_mean_anomaly_thresholds: [f64; 15],
}

impl HyperbolaSolver {
    pub fn new(eccentricity: f64) -> Self {
        let pade_mean_anomaly_thresholds = PADE_ECCENTRIC_ANOMALY_THRESHOLDS.map(|eccentric_anomaly| hyperbolic_kepler_equation(eccentricity, eccentric_anomaly));
        Self { eccentricity, pade_mean_anomaly_thresholds }
    }

    /// Works with all values of mean anomaly
    pub fn solve(&self, mean_anomaly: f64) -> f64 {
        // Solver assumes mean anomaly > 0
        // The equation is symmetric, so for mean anomaly < 0, we just flip the sign o the output
        let ec = self.eccentricity;
        let mh = mean_anomaly.abs();

        let f0 = if mh <= self.pade_mean_anomaly_thresholds[0] {
            // For mh < 5 we use a 'piecewise pade approximation' to get the starting estimate
            let mut i = 0;
            while i < self.pade_mean_anomaly_thresholds.len()-1 && mh < self.pade_mean_anomaly_thresholds[i+1] {
                i += 1;
            }
            let a = PADE_ORDERS[i];
            let coefficients = pade_approximation(ec, mh, a);

            solve_cubic(coefficients, mh, ec) + a

        } else {
            // For mh >= 5, we can use this... thing that I copied from the above paper
            // I have no idea how it works, but it works very very very well
            let fa = f64::ln(2.0 * mh / ec);
            let ca = 0.5 * (2.0 * mh / ec + ec / (2.0 * mh));
            let sa = 0.5 * (2.0 * mh / ec - ec / (2.0 * mh));
            let top = 6.0 * (ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0) 
                + 3.0 * (ec * sa / (ec * ca - 1.0)) * ((ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0)).powi(2);
            let bottom = 6.0 + 6.0 * (ec * sa / (ec * ca - 1.0)) * ((ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0))
                + (ec * ca / (ec * ca - 1.0)) * ((ec.powi(2) / (4.0 * mh) + fa) / (ec * ca - 1.0)).powi(2);
            let delta = top / bottom;
            fa + delta
        };

        // Halley method
        let f = ec * f0.sinh() - f0 - mh;
        let f_prime = ec * f0.cosh() - 1.0;
        let f_prime_prime = f_prime + 1.0;
        let f1 = f0 - (2.0 * f / f_prime) / (2.0 - f * f_prime_prime / f_prime.powi(2));

        // Alternative: Newton-Raphson
        //let f1 = f0 - f / f_prime;

        f1 * mean_anomaly.signum()
    }
}

#[cfg(test)]
mod test {
    use crate::bisection::bisection;

    use super::HyperbolaSolver;

    fn solve_with_bisection(e: f64, m: f64) -> f64 {
        // We don't care about speed here, so just use as wide a range as possible
        let f = |eccentric_anomaly: f64| e * f64::sinh(eccentric_anomaly) - eccentric_anomaly - m;
        bisection(&f, -100000.0, 100000.0)
    }

    #[test]
    fn test_hyperbola() {
        let eccentricites: Vec<f64> = (1..999)
            .map(|x| 1.0 + f64::powi(x as f64, 2) / 1000.0)
            .collect();
        let mean_anomalies: Vec<f64> = (0..10000)
            .map(|x| f64::powi(x as f64, 2) / 10000.0)
            .collect();
            
        for e in &eccentricites {
            let solver = HyperbolaSolver::new(*e);
            for m in &mean_anomalies {
                let expected = solve_with_bisection(*e, *m);
                let actual = solver.solve(*m);
                let difference = if actual.abs() < 1.0e-5 { expected - actual } else { (expected - actual) / actual }.abs();
                if difference > 1.0e-4 {
                    dbg!(expected, actual, difference, e, m);
                    panic!()
                }
            }
        }
    }
}