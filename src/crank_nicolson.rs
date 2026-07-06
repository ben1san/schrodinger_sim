use num_complex::Complex64;

fn solve_tdma(a: Complex64, b: &[Complex64], c: Complex64, d: &[Complex64]) -> Vec<Complex64> {
    let n = d.len();
    let mut c_prime = vec![Complex64::new(0.0, 0.0); n];
    let mut d_prime = vec![Complex64::new(0.0, 0.0); n];
    let mut x = vec![Complex64::new(0.0, 0.0); n];

    c_prime[0] = c / b[0];
    d_prime[0] = d[0] / b[0];

    for i in 1..n {
        let denom = b[i] - a * c_prime[i - 1];
        if i < n - 1 {
            c_prime[i] = c / denom;
        }
        d_prime[i] = (d[i] - a * d_prime[i - 1]) / denom;
    }

    x[n - 1] = d_prime[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    x
}

pub fn evolve_step(psi: &[Complex64], potential: &[f64], dt: f64, dx: f64) -> Vec<Complex64> {
    let n = psi.len();
    let i_comp = Complex64::new(0.0, 1.0);
    let r = i_comp * dt / (4.0 * dx * dx);

    // ディリクレ境界条件では両端(0とn-1)は常に0に固定されるため、
    // 動的に解くべき方程式のサイズは (n-2) となります。
    let inner_n = n - 2;
    let mut d = vec![Complex64::new(0.0, 0.0); inner_n];
    let mut b = vec![Complex64::new(0.0, 0.0); inner_n];

    for i in 1..n - 1 {
        let diag_b = Complex64::new(1.0, 0.0) - i_comp * dt * potential[i] / 2.0 - 2.0 * r;

        // 配列のインデックスを 0 から詰めて格納する (i-1)
        d[i - 1] = r * psi[i - 1] + diag_b * psi[i] + r * psi[i + 1];
        b[i - 1] = Complex64::new(1.0, 0.0) + i_comp * dt * potential[i] / 2.0 + 2.0 * r;
    }

    // サイズ (n-2) の方程式として TDMA を解く
    let inner_psi_next = solve_tdma(-r, &b, -r, &d);

    // 結果を元のサイズ n の配列に書き戻す
    let mut psi_next = vec![Complex64::new(0.0, 0.0); n];
    for i in 1..n - 1 {
        // 両端(0とn-1)は初期値の 0.0+0.0i のまま維持され、内側だけが更新される
        psi_next[i] = inner_psi_next[i - 1];
    }

    psi_next
}

pub fn evolve_step_pbc(psi: &[Complex64], potential: &[f64], dt: f64, dx: f64) -> Vec<Complex64> {
    let n = psi.len();
    let i_comp = Complex64::new(0.0, 1.0);
    let r = i_comp * dt / (4.0 * dx * dx);

    // --- 1. 右辺ベクトル d の計算 (0 と N-1 が隣接するようにラップアラウンド) ---
    let mut d = vec![Complex64::new(0.0, 0.0); n];
    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = if i == n - 1 { 0 } else { i + 1 };

        let diag_b = Complex64::new(1.0, 0.0) - i_comp * dt * potential[i] / 2.0 - 2.0 * r;
        d[i] = r * psi[prev] + diag_b * psi[i] + r * psi[next];
    }

    // --- 2. 左辺の対角成分 b (diag_a) の構築 ---
    let mut b = vec![Complex64::new(0.0, 0.0); n];
    for i in 0..n {
        b[i] = Complex64::new(1.0, 0.0) + i_comp * dt * potential[i] / 2.0 + 2.0 * r;
    }

    // --- 3. Sherman-Morrisonの公式による巡回三重対角行列の求解 ---
    let a = -r; // 下側対角成分
    let c = -r; // 上側対角成分

    // A' の対角成分を修正 (角の要素を分離)
    let mut b_prime = b.clone();
    b_prime[0] = b_prime[0] - a;
    b_prime[n - 1] = b_prime[n - 1] - c;

    // u ベクトルの作成
    let mut u = vec![Complex64::new(0.0, 0.0); n];
    u[0] = a;
    u[n - 1] = c;

    // A'y = d と A'z = u をそれぞれ標準TDMAで解く
    let y = solve_tdma(a, &b_prime, c, &d);
    let z = solve_tdma(a, &b_prime, c, &u);

    // v = [1, 0, ..., 0, 1]^T との内積計算
    let v_dot_y = y[0] + y[n - 1];
    let v_dot_z = z[0] + z[n - 1];

    // 解の合成
    let factor = v_dot_y / (Complex64::new(1.0, 0.0) + v_dot_z);
    let mut psi_next = vec![Complex64::new(0.0, 0.0); n];
    for i in 0..n {
        psi_next[i] = y[i] - factor * z[i];
    }

    psi_next
}
