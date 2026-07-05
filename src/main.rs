use num_complex::Complex64;

fn main() {
    const DT: f64 = 0.1;
    const DX: f64 = 0.1;
    const N: usize = 1000;
    const TIMESTEPS: usize = 1000;
    let potential = vec![0.0; N];

    let mut psi = vec![Complex64::new(0.0, 0.0); N];
    let i_comp = Complex64::new(0.0, 1.0);

    let x0 = (N as f64 * DX) / 4.0;
    let sigma = 5.0;
    let k0 = 2.0;

    for i in 1..N - 1 {
        let x = i as f64 * DX;
        let envelope = f64::exp(-0.5 * ((x - x0) / sigma).powi(2));
        let phase = k0 * x;
        psi[i] = Complex64::new(envelope * phase.cos(), envelope * phase.sin())
    }

    let r = i_comp * DT / (4.0 * DX * DX);

    for _step in 0..TIMESTEPS {
        let mut d = vec![Complex64::new(0.0, 0.0); N];
        for i in 1..N - 1 {
            let diag_b = Complex64::new(1.0, 0.0) - i_comp * DT * potential[i] / 2.0 - 2.0 * r;
            d[i] = r * psi[i - 1] + diag_b * psi[i] + r * psi[i + 1]
        }

        let mut c_prime = vec![Complex64::new(0.0, 0.0); N];
        let mut d_prime = vec![Complex64::new(0.0, 0.0); N];
        let mut psi_next = vec![Complex64::new(0.0, 0.0); N];

        for i in 1..N - 1 {
            let diag_a = Complex64::new(1.0, 0.0) + i_comp * DT * potential[i] / 2.0 + 2.0 * r;
            let a = -r;
            let c = -r;

            if i == 1 {
                c_prime[i] = c / diag_a;
                d_prime[i] = d[i] / diag_a;
            } else {
                c_prime[i] = c / (diag_a - c_prime[i - 1] * a);
                d_prime[i] = (d[i] - d_prime[i - 1] * a) / (diag_a - c_prime[i - 1] * a);
            }
        }

        for i in (1..N - 1).rev() {
            if i == N - 2 {
                psi_next[i] = d_prime[i];
            } else {
                psi_next[i] = d[i] - c_prime[i] * psi_next[i + 1]
            }
        }

        psi = psi_next;
    }
}
