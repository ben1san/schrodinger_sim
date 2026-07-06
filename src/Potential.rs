/// 自由空間 (V = 0) のポテンシャルを生成
pub fn free_space(n: usize) -> Vec<f64> {
    vec![0.0; n]
}

/// ポテンシャル障壁を生成（トンネル効果のシミュレーション用）
/// start_x: 障壁の開始位置, width: 障壁の幅, height: 障壁の高さ(エネルギー)
pub fn barrier(n: usize, dx: f64, start_x: f64, width: f64, height: f64) -> Vec<f64> {
    let mut v = vec![0.0; n];
    for i in 0..n {
        let x = i as f64 * dx;
        if x >= start_x && x <= start_x + width {
            v[i] = height;
        }
    }
    v
}
