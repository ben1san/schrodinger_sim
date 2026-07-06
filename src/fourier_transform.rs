use num_complex::Complex64;
use std::f64::consts::PI;

/// 位置空間の波動関数 psi(x) を離散フーリエ変換し、
/// 運動量空間の波動関数 phi(p) を計算して返す
pub fn momentum_distribution(psi: &[Complex64], dx: f64) -> Vec<Complex64> {
    let n = psi.len();
    let mut phi = vec![Complex64::new(0.0, 0.0); n];

    // 運動量空間の刻み幅 (フーリエ変換の性質から自動的に決まる)
    let dp = 2.0 * PI / (n as f64 * dx);

    for k in 0..n {
        // インデックス k=0~N を、負から正の運動量 p にマッピングする
        // これにより、配列の中央(N/2)が運動量0になります。
        let p = dp * (k as f64 - (n as f64) / 2.0);

        let mut phi_k = Complex64::new(0.0, 0.0);

        for j in 0..n {
            let x = j as f64 * dx;
            // e^{-i p x}
            let phase = Complex64::new(0.0, -p * x);

            // 積分(総和)の実行: psi(x) * e^{-i p x} * dx
            phi_k += psi[j] * phase.exp() * dx;
        }

        // 規格化: 1/sqrt(2π)
        phi[k] = phi_k / (2.0 * PI).sqrt();
    }

    phi
}
