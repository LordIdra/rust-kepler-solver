const BISECTION_ITERATIONS: usize = 64;

pub fn bisection(function: &impl Fn(f64) -> f64, min: f64, max: f64) -> f64 {
    let mut low = min;
    let mut high = max;
    let mut mid = (min + max) / 2.0;
    for _ in 0..BISECTION_ITERATIONS {
        if function(mid).is_sign_positive() && function(low).is_sign_positive() || function(mid).is_sign_negative() && function(low).is_sign_negative() {
            low = mid;
        } else {
            high = mid;
        }
        mid = (low + high) / 2.0;
    }
    mid
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use crate::bisection::bisection;

    #[test]
    fn test_bisection_1() {
        let function = |x: f64| f64::sin(x);
        assert!(bisection(&function, -PI/2.0, PI/2.0).abs() < 1.0e-2);
    }

    #[test]
    fn test_bisection_2() {
        let function = |x: f64| x.powi(2) - 4.0;
        assert!(bisection(&function, 0.0,  10.0).abs() - 2.0 < 1.0e-2);
        assert!(bisection(&function, -10.0, 0.0).abs() - 2.0 < 1.0e-2);
    }
}