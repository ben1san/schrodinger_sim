use num_complex::Complex64;

pub fn evolve_step(psi: &[Complex64], potential: &[f64], dt: f64, dx: f64) -> Vec<Complex64> {
    let n = psi.len();
    let i_comp = Complex64::new(0.0, 1.0);
    let r = i_comp * dt / (4.0 * dx * dx);

    let mut d = vec![Complex64::new(0.0, 0.0); n];
    for i in 1..n - 1 {
        let diag_b = Complex64::new(1.0, 0.0) - i_comp * dt * potential[i] / 2.0 - 2.0 * r;
        d[i] = r * psi[i - 1] + diag_b * psi[i] + r * psi[i + 1];
    }

    let mut c_prime = vec![Complex64::new(0.0, 0.0); n];
    let mut d_prime = vec![Complex64::new(0.0, 0.0); n];
    let mut psi_next = vec![Complex64::new(0.0, 0.0); n];

    for i in 1..n - 1 {
        let diag_a = Complex64::new(1.0, 0.0) + i_comp * dt * potential[i] / 2.0 + 2.0 * r;
        let a = -r;
        let c = -r;

        if i == 1 {
            c_prime[i] = c / diag_a;
            d_prime[i] = d[i] / diag_a;
        } else {
            let denom = diag_a - a * c_prime[i - 1];
            c_prime[i] = c / denom;
            d_prime[i] = (d[i] - a * d_prime[i - 1]) / denom;
        }
    }

    for i in (1..n - 1).rev() {
        if i == n - 2 {
            psi_next[i] = d_prime[i];
        } else {
            psi_next[i] = d_prime[i] - c_prime[i] * psi_next[i + 1];
        }
    }

    psi_next
}
