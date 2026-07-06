mod crank_nicolson;
mod fourier_transform;
mod initial_state;
mod potential;

use std::fs::File;
use std::io::{BufWriter, Write};
fn main() {
    const DT: f64 = 0.01;
    const DX: f64 = 0.1;
    const N: usize = 1000;
    const TIME_STEPS: usize = 3000;

    let potential = potential::barrier(N, DX, 45.0, 55.0, 2.0);
    let mut psi = initial_state::normal_distr(N, DX);

    let mut writer_position = BufWriter::new(File::create("data/probability_density.csv").unwrap());
    let mut writer_momentum = BufWriter::new(File::create("data/momentum_density.csv").unwrap());

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
            writer_position
                .write_all(line.as_bytes())
                .expect("書き込みエラー");

            let phi = fourier_transform::momentum_distribution(&psi, DX);
            let mut line = String::new();
            for i in 0..N {
                let density = phi[i].norm_sqr();

                line.push_str(&format!("{:.6}", density));

                if i < N - 1 {
                    line.push(',');
                }
            }
            line.push('\n');
            writer_momentum
                .write_all(line.as_bytes())
                .expect("書き込みエラー");
        }
    }
}
