use pa_types::cigar::*;

pub fn run() {
    println!("Hello, world!");
}

#[inline(always)]
fn collect_match_positions(cigar: &Cigar) -> Vec<usize> {
    let mut match_positions = Vec::new();
    let mut current_pos: usize = 0;
    for elem in &cigar.ops {
        let count = elem.cnt as usize;
        if elem.op == CigarOp::Match {
            for offset in 0..count {
                match_positions.push(current_pos + offset);
            }
        }
        current_pos += count;
    }
    match_positions
}

#[inline(always)]
fn compute(cigar: &Cigar, k: usize, lambda_decay: f64) -> f64 {
    assert!(k > 0, "k must be greater than 0");

    let match_pos = collect_match_positions(cigar);

    let m = match_pos.len();
    if m < k {
        return 0.0;
    }

    // Base for subsequences of length 1: λ^1 per match position
    let mut dp_prev = vec![lambda_decay; m];
    if k == 1 {
        return dp_prev.iter().sum::<f64>();
    }

    // Precompute λ^pos[0] and jump multipliers λ^{pos[i+1]-pos[i]}
    let lam_start = lambda_decay.powi(match_pos[0] as i32);
    let mut lam_jumps: Vec<f64> = Vec::with_capacity(m.saturating_sub(1));
    for i in 0..(m.saturating_sub(1)) {
        let delta = (match_pos[i + 1] - match_pos[i]) as i32;
        lam_jumps.push(lambda_decay.powi(delta));
    }

    // dp_l[i] = λ^{pos[i]} * sum_{j < i} (dp_{l-1}[j] * λ^{-pos[j]})
    let mut dp = vec![0.0; m];
    for _ in 2..=k {
        let mut running = 0.0;
        let mut lam_cur = lam_start;
        let mut lam_cur_inv = 1.0 / lam_cur;

        for i in 0..m {
            dp[i] = lam_cur * running;
            running += dp_prev[i] * lam_cur_inv;

            if i + 1 < m {
                let jump = lam_jumps[i];
                lam_cur *= jump;
                lam_cur_inv /= jump;
            }
        }
        dp_prev.copy_from_slice(&dp);
    }

    dp_prev.iter().sum::<f64>()
}

pub struct Lodhi {
    k: usize,
    lambda_decay: f64,
    match_pos: Vec<usize>,
    lam_jumps: Vec<f64>,
    dp_prev: Vec<f64>,
    dp: Vec<f64>,
}

impl Lodhi {
    pub fn new(k: usize, lambda_decay: f64) -> Self {
        Self {
            k,
            lambda_decay,
            match_pos: Vec::new(),
            lam_jumps: Vec::new(),
            dp_prev: Vec::new(),
            dp: Vec::new(),
        }
    }

    #[inline(always)]
    pub fn compute(&mut self, cigar: &Cigar) -> f64 {
        assert!(self.k > 0, "k must be greater than 0");

        self.match_pos.clear();
        let mut current_pos: usize = 0;
        for elem in &cigar.ops {
            let count = elem.cnt as usize;
            if elem.op == CigarOp::Match {
                for offset in 0..count {
                    self.match_pos.push(current_pos + offset);
                }
            }
            current_pos += count;
        }

        let m = self.match_pos.len();
        if m < self.k {
            return 0.0;
        }

        // Prepare buffers
        self.dp_prev.resize(m, 0.0);
        self.dp.resize(m, 0.0);
        self.lam_jumps.resize(m.saturating_sub(1), 0.0);

        // Base case: λ per match position
        self.dp_prev.fill(self.lambda_decay);
        if self.k == 1 {
            return self.dp_prev.iter().sum::<f64>();
        }

        // Precompute λ^pos[0] and jump multipliers λ^{pos[i+1]-pos[i]}
        let lam_start = self.lambda_decay.powi(self.match_pos[0] as i32);
        for i in 0..(m.saturating_sub(1)) {
            let delta = (self.match_pos[i + 1] - self.match_pos[i]) as i32;
            self.lam_jumps[i] = self.lambda_decay.powi(delta);
        }

        for _ in 2..=self.k {
            let mut running = 0.0;
            let mut lam_cur = lam_start;
            let mut lam_cur_inv = 1.0 / lam_cur;

            for i in 0..m {
                self.dp[i] = lam_cur * running;
                running += self.dp_prev[i] * lam_cur_inv;

                if i + 1 < m {
                    let jump = self.lam_jumps[i];
                    lam_cur *= jump;
                    lam_cur_inv /= jump;
                }
            }
            self.dp_prev.copy_from_slice(&self.dp);
        }

        self.dp_prev.iter().sum::<f64>()
    }
}

mod tests {
    use super::*;
    use pa_types::*;
    use std::time::Instant;

    #[test]
    fn test_compute() {
        //2=1D3=
        let ops: Vec<CigarOp> = vec![
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Del,
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Match,
        ];

        let cigar = Cigar::from_ops(ops.into_iter());
        let start_time = Instant::now();
        let result = compute(&cigar, 3, 0.5);
        let end_time = Instant::now();
        println!("result: {result}");
        println!("time: {:?}", end_time.duration_since(start_time));
        assert_eq!(result, 0.421875);
    }

    #[test]
    fn test_lodhi_compute_reuse() {
        let ops: Vec<CigarOp> = vec![
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Del,
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Match,
        ];
        let cigar = Cigar::from_ops(ops.into_iter());

        let mut lodhi = Lodhi::new(3, 0.5);
        let r3 = lodhi.compute(&cigar);
        assert_eq!(r3, 0.421875);

        // Reuse allocations across a different cigar
        let ops2: Vec<CigarOp> = vec![
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Del,
            CigarOp::Match,
        ];
        let cigar2 = Cigar::from_ops(ops2.into_iter());
        let start_time = Instant::now();
        let r3_b = lodhi.compute(&cigar2);
        let end_time = Instant::now();
        println!(
            "Time re-using alloc: {:?}",
            end_time.duration_since(start_time)
        );
        assert_eq!(r3_b, compute(&cigar2, 3, 0.5));
    }

    #[test]
    fn compute_paper_cigar() {
        let mut lodhi = Lodhi::new(3, 0.5);

        let ops1: Vec<CigarOp> = vec![
            CigarOp::Sub,
            CigarOp::Sub,
            CigarOp::Match,
            CigarOp::Match,
            CigarOp::Match,
        ];
        let cigar1 = Cigar::from_ops(ops1.into_iter());
        let start_time = Instant::now();
        let s1 = lodhi.compute(&cigar1);
        let end_time = Instant::now();
        println!("s1: {s1}");
        println!("time: {:?}", end_time.duration_since(start_time));

        let ops2: Vec<CigarOp> = vec![
            CigarOp::Match,
            CigarOp::Sub,
            CigarOp::Match,
            CigarOp::Sub,
            CigarOp::Match,
        ];
        let cigar2 = Cigar::from_ops(ops2.into_iter());
        let start_time = Instant::now();
        let s2 = lodhi.compute(&cigar2);
        let end_time = Instant::now();
        println!("s2: {s2}");
        println!("time: {:?}", end_time.duration_since(start_time));
    }
}
