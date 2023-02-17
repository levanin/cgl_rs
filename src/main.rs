use std::env;
use ark_ff::Field;
use ark_ff::fields::{Fp320, MontBackend, Fp2Config, Fp2, MontFp};
// init base field
#[derive(ark_ff::MontConfig)]
#[modulus = "5210644015679228794060694325390955853397149309953825381775591280356090833797119"]
#[generator = "17"]
pub struct FqConfig;
pub type Fq = Fp320<MontBackend<FqConfig, 5>>;

pub const FQ_ONE: Fq = MontFp!("1");
pub const FQ_ZERO: Fq = MontFp!("0");
pub const FQ_MINUS_ONE: Fq = MontFp!("-1");
// init quadratic extension field
pub struct Fq2Config;
pub type Fq2 = Fp2<Fq2Config>;

impl Fp2Config for Fq2Config{
    type Fp = Fq;
    #[rustfmt::skip]
    const NONRESIDUE: Fq = MontFp!("-1");
    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[
        // Fq(-1)**(((q^0) - 1) / 2)
        MontFp!("1"),
        // Fq(-1)**(((q^1) - 1) / 2)
        MontFp!("-1"),
    ];
    #[inline(always)]
    fn mul_fp_by_nonresidue_in_place(fe: &mut Self::Fp) -> &mut Self::Fp {
        *fe *= Self::NONRESIDUE;
        fe
    }

    /// A specializable method for setting `y = x + NONRESIDUE * y`.
    /// This allows for optimizations when the non-residue is
    /// canonically negative in the field.
    #[inline(always)]
    fn mul_fp_by_nonresidue_and_add(y: &mut Self::Fp, x: &Self::Fp) {
        Self::mul_fp_by_nonresidue_in_place(y);
        *y += x;
    }

    /// A specializable method for computing x + mul_fp_by_nonresidue(y) + y
    /// This allows for optimizations when the non-residue is not -1.
    #[inline(always)]
    fn mul_fp_by_nonresidue_plus_one_and_add(y: &mut Self::Fp, x: &Self::Fp) {
        let old_y = *y;
        Self::mul_fp_by_nonresidue_and_add(y, x);
        *y += old_y;
    }

    /// A specializable method for computing x - mul_fp_by_nonresidue(y)
    /// This allows for optimizations when the non-residue is
    /// canonically negative in the field.
    #[inline(always)]
    fn sub_and_mul_fp_by_nonresidue(y: &mut Self::Fp, x: &Self::Fp) {
        *y = *x - Self::mul_fp_by_nonresidue_in_place(y);
    }
}
pub const FQ2_ZERO: Fq2 = Fq2::new(FQ_ZERO, FQ_ZERO);
pub const FQ2_ONE: Fq2 = Fq2::new(FQ_ONE, FQ_ZERO);
pub const FQ2_MINUS_ONE: Fq2 = Fq2::new(FQ_MINUS_ONE, FQ_ZERO);
pub const C0: Fq2 = Fq2::new(MontFp!("1488"), MontFp!("0"));
pub const C1: Fq2 = Fq2::new(MontFp!("-162000"), MontFp!("0"));
pub const C2: Fq2 = Fq2::new(MontFp!("40773375"), MontFp!("0"));
pub const C3: Fq2 = Fq2::new(MontFp!("8748000000"), MontFp!("0"));
// Idk why but I can't figure out how to define a const which is an inverse of something. Below is 2^-1 in Fq2
pub const C4: Fq2 = Fq2::new(MontFp!("2605322007839614397030347162695477926698574654976912690887795640178045416898560"), MontFp!("0"));
pub const C5: Fq2 = Fq2::new(MontFp!("4"), MontFp!("0"));

fn walk_step(j_i: &Fq2, j_i_minus_one: &Fq2, sign: Fq2) -> (Fq2,Fq2){
    let j_i_sqr: Fq2 = (*j_i).square();
    let j_i_minus_one_sqr: Fq2 = (*j_i_minus_one).square();

    let a_i: Fq2 = -j_i_sqr + C0 * *j_i + C1;
    let b_i: Fq2 = C0*j_i_sqr + C2* *j_i + C3;
    let d_i: Fq2 = (a_i + j_i_minus_one).square() - C5*(b_i + a_i* *j_i_minus_one + j_i_minus_one_sqr);
    let s_i: Fq2 = d_i.sqrt().unwrap();

    (C4*(-a_i - *j_i_minus_one + sign*s_i), *j_i)
}

fn main() {
    let args:Vec<String> = env::args().collect();
    let input:&str = &args[1];
    let mut j_i: Fq2 = Fq2::new(MontFp!("1728"), MontFp!("0"));
    let mut j_i_minus_one: Fq2 = Fq2::new(MontFp!("1728"), MontFp!("0"));
    let mut sign: Fq2;
    for c in input.chars(){
        assert!(c == '0' || c == '1');
        sign = if c == '0'{FQ2_MINUS_ONE} else {FQ2_ONE};
        (j_i, j_i_minus_one) = walk_step(&j_i, &j_i_minus_one, sign);
    }
    if &j_i.c1.to_string() == ""
    {
        println!("{}", &j_i.c0.to_string());
    }
    else
    {
        println!("{} + {}*i", &j_i.c0.to_string(), &j_i.c1.to_string());
    }




}
// prime chosen is 2^256 * 45 - 1