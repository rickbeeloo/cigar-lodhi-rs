use cigar_lodhi_rs::*;
use pa_types::*;

fn main() {
    let mut lodhi = Lodhi::new(3, 0.5);
    let cigar = Cigar::from_ops(
        vec![
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Del,
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Match,
        ]
        .into_iter(),
    );
    let score = lodhi.compute(&cigar);
    println!("Score: {}", score);
}
