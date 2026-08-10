pub fn sinh_cosh(a: f64) -> (f64, f64) {
    let ex = f64::exp(a);
    let enx = f64::exp(-a);
    let sinh = (ex-enx) / 2.0;
    let cosh = (ex+enx) / 2.0;
    (sinh, cosh)
}

pub fn arsinh(x: f64) -> f64 {
    f64::ln(x + f64::sqrt(x.powi(2) + 1.0))
}
