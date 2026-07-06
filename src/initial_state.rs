use num_complex::Complex64;
pub fn normal_distr(n: usize, dx: f64) -> Vec<Complex64> {
    let x0 = (n as f64 * dx) / 4.0;
    let sigma = 5.0;
    let k0 = 2.0;
    let mut psi_init = vec![Complex64::new(0.0, 0.0); n];

    for i in 1..n - 1 {
        let x = i as f64 * dx;
        let envelope = f64::exp(-0.5 * ((x - x0) / sigma).powi(2));
        let phase = k0 * x;
        psi_init[i] = Complex64::new(envelope * phase.cos(), envelope * phase.sin());
    }

    psi_init
}
