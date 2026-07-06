mod crank_nicolson;
mod initial_state;
mod potential;

use std::fs::File;
use std::io::{BufWriter, Write};
fn main() {
    const DT: f64 = 0.01;
    const DX: f64 = 0.1;
    const N: usize = 1000;
    const TIME_STEPS: usize = 3000;

    let potential = potential::free_space(N);
    let mut psi = initial_state::normal_distr(N, DX);

    let filename = if std::path::Path::new("../data").exists() {
        "../data/probability_density.csv"
    } else {
        "data/probability_density.csv"
    };
    let file = File::create(filename).expect("ファイルの作成に失敗しました");
    let mut writer = BufWriter::new(file);

    for step in 0..TIME_STEPS {
        psi = crank_nicolson::evolve_step(&psi, &potential, DT, DX);

        if step % 10 == 0 {
            let mut line = String::new();
            for i in 0..N {
                let density = psi[i].norm_sqr();

                line.push_str(&format!("{:.6}", density));

                if i < N - 1 {
                    line.push(',');
                }
            }
            line.push('\n');
            writer.write_all(line.as_bytes()).expect("書き込みエラー");
        }
    }
}
