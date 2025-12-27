//! Implementation of the SIP standard.

use crate::{CustomFloat, ImgXY};
use std::ops::RangeInclusive;

/// SIP Polynomial coefficient.
/// In the polynomial, coefficient must be ordered like this:
/// * `0_0, 0_1, 0_2, 0_3, 1_0, 1_1, 1_2, 2_0, 2_1, 3_0`
/// * in which `p_q` correspond to the polynomial part `coeff_p_q * u^p * v^q`
/// * Given an order `n`, the size of the array must be `n(n+1)/2`.
#[derive(Debug, Clone)]
pub struct SipCoeff {
    /// Computed order of the polynomial
    order: u16,

    /// Polynomials coefficient matrix
    c: Box<[f64]>,
}

impl SipCoeff {
    /// # Param
    /// * `c`: array polynomial coefficients of size `n(n+1)/2`
    #[must_use]
    pub fn new(c: Box<[f64]>) -> Self {
        let n_coeff = c.len();

        let order = ((((n_coeff * 8) + 1) as f64).sqrt() as u16 - 1) / 2;

        debug_assert_eq!(order * (order + 1) / 2, (c.len() as u16));
        Self { order, c }
    }

    /// Returns the value of the polynomial, evaluated in `(u, v)`.
    #[must_use]
    pub fn p(&self, u: f64, v: f64) -> f64 {
        let mut k = 0;
        let mut p = 0_f64;
        let mut v_pow = 1.0;
        let mut u_pow: f64;
        // loop over v^i
        for i in 0..self.order {
            u_pow = 1.0;
            // loop over u^j
            for _ in 0..self.order - i {
                p += u_pow * v_pow * unsafe { self.c.get_unchecked(k) };
                k += 1;
                u_pow *= u;
            }
            v_pow *= v;
        }
        debug_assert_eq!(k, self.c.len(), "Should have iterated over all of c");
        p
    }

    /// Returns the value of the `dp/du`, evaluated in `(u, v)`.
    #[must_use]
    pub fn dpdu(&self, u: f64, v: f64) -> f64 {
        let mut k = 0;
        let mut p = 0_f64;
        let mut v_pow = 1.0;
        let mut u_pow: f64;
        for i in 0..self.order {
            u_pow = 1.0;
            k += 1;
            for j in 1..self.order - i {
                p += u_pow * v_pow * unsafe { self.c.get_unchecked(k) } * f64::from(j);
                k += 1;
                u_pow *= u;
            }
            v_pow *= v;
        }
        p
    }

    /// Returns the value of the `dp/dv`, evaluated in `(u, v)`.
    #[must_use]
    pub fn dpdv(&self, u: f64, v: f64) -> f64 {
        let mut k = (self.order) as usize;
        let mut p = 0_f64;
        let mut v_pow: f64 = 1.0;
        let mut u_pow: f64;
        for i in 1..self.order {
            u_pow = 1.0;
            for _ in 0..self.order - i {
                p += u_pow * v_pow * unsafe { self.c.get_unchecked(k) } * f64::from(i);
                k += 1;
                u_pow *= u;
            }
            v_pow *= v;
        }
        p
    }
}

/// SIP (un)projection coefficients for 1st and 2nd axis
#[derive(Debug, Clone)]
pub struct SipAB {
    /// Polynomials coefficient matrix on the 1st axis.
    a: SipCoeff,

    /// Polynomials coefficient matrix on the2ndt axis.
    b: SipCoeff,
}

impl SipAB {
    /// # Params
    /// * `a`: 1st axis SIP coefficients
    /// * `b`: 2nd axis SIP coefficients
    #[must_use]
    pub fn new(a: SipCoeff, b: SipCoeff) -> Self {
        Self { a, b }
    }
}

/// For the SIP convention, see
/// "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
/// in the proceedings of ADASS XIV (2005).
#[derive(Debug, Clone)]
pub struct Sip {
    /// Projection coefficient.
    ab_proj: SipAB,

    /// De-projection coefficients (if any).
    ab_deproj: Option<SipAB>,

    /// Approximate bounds of the 1st axis domain of validity.
    /// * `start`: `-(CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
    /// * `end`: `(NAXIS1 - CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
    u: RangeInclusive<f64>,

    /// Approximate bounds of the 2nd axis domain of validity.
    /// * `start`: `-(CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
    /// * `end`: `(NAXIS2 - CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
    v: RangeInclusive<f64>,

    fuv: RangeInclusive<f64>,
    guv: RangeInclusive<f64>,

    /// Number of iteration of the mutli-variate Newton-Raphson method (if no unproj polynomial).
    n_iter: u8,

    /// Precision used in the mutli-variate Newton-Raphson method (if no unproj polynomial).
    eps: f64,
}

impl Sip {
    /// Implements the SIP convention with the given polynomial coefficients.
    /// # Params
    /// * `ab_proj`: SIP coefficients for the projection on the 1st and 2nd axis
    /// * `ab_deproj`: SIP coefficients for the deprojection on the 1st and 2nd axis (if any)
    /// * `u`: 1st axis domain of validity, e.g. `[-CRPIX1..NAXIS1 - CRPIX1]`
    /// * `v`: 2nd axis domain of validity, e.g. `[-CRPIX2..NAXIS2 - CRPIX2]`
    #[must_use]
    pub fn new(
        ab_proj: SipAB,
        ab_deproj: Option<SipAB>,
        u: RangeInclusive<f64>,
        v: RangeInclusive<f64>,
    ) -> Self {
        let t = ab_proj
            .a
            .p(*u.start(), *v.start())
            .min(ab_proj.a.p(*u.start(), *v.end()));
        let fuv_min = ab_proj.a.p(*u.start(), 0.0).min(t);
        let t = ab_proj
            .a
            .p(*u.end(), *v.start())
            .max(ab_proj.a.p(*u.end(), *v.end()));
        let fuv_max = ab_proj.a.p(*u.end(), 0.0).max(t);

        let t = ab_proj
            .b
            .p(*u.start(), *v.start())
            .min(ab_proj.b.p(*u.end(), *v.start()));
        let guv_min = ab_proj.b.p(0.0, *v.start()).min(t);
        let t = ab_proj
            .b
            .p(*u.start(), *v.end())
            .max(ab_proj.b.p(*u.end(), *v.end()));
        let guv_max = ab_proj.b.p(0.0, *v.end()).max(t);

        Self {
            ab_proj,
            ab_deproj,
            u,
            v,
            fuv: fuv_min..=fuv_max,
            guv: guv_min..=guv_max,
            n_iter: 20,
            eps: 1.0e-9,
        }
    }

    /// Does SIP contain polynomial deprojection
    #[must_use]
    pub fn has_polynomial_deproj(&self) -> bool {
        self.ab_deproj.is_some()
    }

    #[must_use]
    pub(crate) fn f(&self, u: f64, v: f64) -> f64 {
        self.ab_proj.a.p(u, v)
    }

    #[must_use]
    pub(crate) fn g(&self, u: f64, v: f64) -> f64 {
        self.ab_proj.b.p(u, v)
    }

    #[must_use]
    pub(crate) fn dfdu(&self, u: f64, v: f64) -> f64 {
        self.ab_proj.a.dpdu(u, v)
    }

    #[must_use]
    pub(crate) fn dfdv(&self, u: f64, v: f64) -> f64 {
        self.ab_proj.a.dpdv(u, v)
    }

    #[must_use]
    pub(crate) fn dgdu(&self, u: f64, v: f64) -> f64 {
        self.ab_proj.b.dpdu(u, v)
    }

    #[must_use]
    pub(crate) fn dgdv(&self, u: f64, v: f64) -> f64 {
        self.ab_proj.b.dpdv(u, v)
    }

    #[must_use]
    pub(crate) fn u(&self, fuv: f64, guv: f64) -> Option<f64> {
        self.ab_deproj.as_ref().map(|ab| ab.a.p(fuv, guv))
    }

    #[must_use]
    pub(crate) fn v(&self, fuv: f64, guv: f64) -> Option<f64> {
        self.ab_deproj.as_ref().map(|ab| ab.b.p(fuv, guv))
    }

    #[must_use]
    pub(crate) fn inverse(&self, fuv: f64, guv: f64) -> Option<ImgXY> {
        // uv
        if self.has_polynomial_deproj() {
            let u = self.u(fuv, guv).unwrap();
            let v = self.v(fuv, guv).unwrap();
            Some(ImgXY::new(u, v))
        } else {
            // TODO: Make a grid and a 2-d tree to find the starting point (then multi-variate Newton)
            None
        }
    }

    /// Multi-variate Newton-Raphson:
    /// f1(x1, ..., xn) = 0es006500
    /// ...
    /// fn(x1, ..., xn) = 0
    ///  
    /// x = (x1, ..., xn)
    /// f = (f1(x), ... fn(x))
    ///  
    /// x = x - J^-1 f
    ///
    ///  With J = df1/dx1 ... df1/dxn
    /// ...     ... ...
    /// dfn/dx1 ... dfn/dxn
    ///
    /// 2d case: M = ab => M^-1 = 1/(ad-bc)  d -b
    /// cd                     -c  a
    ///
    #[must_use]
    pub fn bivariate_newton(&self, fuv: f64, guv: f64) -> Option<ImgXY> {
        // Check input values are in the domain of validity
        if self.fuv.contains(&fuv) && self.guv.contains(&guv) {
            // Initial guess
            let mut u = fuv;
            let mut v = guv;
            // Initial values
            let mut f = self.f(u, v) - fuv;
            let mut g = self.g(u, v) - guv;
            // Bivariate Newton's method
            let eps2 = self.eps.pow2();
            let mut norm2 = f.pow2() + g.pow2();
            let mut i = 0;
            while i < self.n_iter && norm2 < eps2 {
                let a = self.dfdu(u, v);
                let b = self.dfdv(u, v);
                let c = self.dgdu(u, v);
                let d = self.dgdv(u, v);
                let det = 1.0 / (a * d - b * c);
                u -= det * (f * d - g * b);
                v -= det * (g * a - f * c);
                f = self.f(u, v) - fuv;
                g = self.g(u, v) - guv;
                norm2 = f.pow2() + g.pow2();
                i += 1;
            }
            // Check that the result is in the domain of validity
            if self.u.contains(&u) && self.v.contains(&v) {
                // All good :)
                Some(ImgXY::new(u, v))
            } else {
                // TODO: look for a different initial guess by making a grid of 300x300 and put results
                // in a kd-tree. Then look at the nearest neighbour: it is the initial guess.
                // And redo Newton.
                None
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sip::{Sip, SipAB, SipCoeff};

    #[test]
    fn test_sip() {
        // taken from table 1 of
        // "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
        // in the proceedings of ADASS XIV (2005).
        //
        // CTYPE1 ’RA---TAN-SIP’
        // CTYPE2 ’DEC--TAN-SIP’
        // CRPIX1 2048.0
        // CRPIX2 1024.0
        // CRVAL1 5.6260667398471
        // CRVAL2 -72.076963036772
        // CD1_1 -7.8481866550866E-06
        // CD2_1 1.1406694624771E-05
        // CD1_2 1.0939720432379E-05
        // CD2_2 8.6942510845452E-06

        // A_0_0 0.0
        // A_0_1 0.0
        // A_0_2 2.1634068532689E-06
        // A_0_3 1.0622437604068E-11
        // A_0_4 1.4075878614807E-14
        // A_1_0 0.0
        // A_1_1 -5.194753640575E-06
        // A_1_2 -5.2797808038221E-10
        // A_1_3 -1.9317154005522E-14
        // A_2_0 8.543473309812E-06
        // A_2_1 -4.4012735467525E-11
        // A_2_2 3.767898933666E-14
        // A_3_0 -4.7518233007536E-10
        // A_3_1 5.0860953083043E-15
        // A_4_0 2.5776347115304E-14

        // B_0_0 0.0
        // B_0_1 0.0
        // B_0_2 -7.2299995118730E-06
        // B_0_3 -4.2102920235938E-10
        // B_0_4 6.5531313110898E-16
        // B_1_0 0.0
        // B_1_1 6.1778338717084E-06
        // B_1_2 -6.7603466821178E-11
        // B_1_3 1.3892905568706E-14
        // B_2_0 -1.7442694174934E-06
        // B_2_1 -5.1333879897858E-10
        // B_2_2 -2.9648166208490E-14
        // B_3_0 8.5722142612681E-11
        // B_3_1 -2.0749495718513E-15
        // B_4_0 -1.812610418272E-14

        // A_ORDER 4
        // B_ORDER 4

        // In the polynomial, coefficient must be ordered like this:
        // * `0_0, 0_1, 0_2, 0_3, 1_0, 1_1, 1_2, 2_0, 2_1, 3_0`
        let a_coef = [
            0.0,
            0.0,
            2.1634068532689E-06,
            1.0622437604068E-11,
            1.4075878614807E-14,
            0.0,
            -5.194753640575E-06,
            -5.2797808038221E-10,
            -1.9317154005522E-14,
            8.543473309812E-06,
            -4.4012735467525E-11,
            3.767898933666E-14,
            -4.7518233007536E-10,
            5.0860953083043E-15,
            2.5776347115304E-14,
        ];

        let a_coef = SipCoeff::new(a_coef.into());

        let b_coef = [
            0.0,
            0.0,
            -7.2299995118730E-06,
            -4.2102920235938E-10,
            6.5531313110898E-16,
            0.0,
            6.1778338717084E-06,
            -6.7603466821178E-11,
            1.3892905568706E-14,
            -1.7442694174934E-06,
            -5.1333879897858E-10,
            -2.9648166208490E-14,
            8.5722142612681E-11,
            -2.0749495718513E-15,
            -1.812610418272E-14,
        ];

        let b_coef = SipCoeff::new(b_coef.into());
        let ab_coef = SipAB::new(a_coef, b_coef);

        let sip = Sip::new(ab_coef, None, -2048.0..=2048.0, -1024.0..=1024.0);

        let tmp = sip.f(3.0, 3.0);
        let tmp2 = sip.g(3.0, 3.0);
        let dfdu = sip.dfdu(3.0, 3.0);
        let dfdv = sip.dfdv(3.0, 3.0);

        let da = sip.f(3.0, 3.0);
        let db = sip.f(3.0 + 1e-9, 3.0);
        let dc = sip.f(3.0, 3.0 + 1e-9);
        let dfdu_approx = (db - da) / 1e-9;
        let dfdv_approx = (dc - da) / 1e-9;

        dbg!(tmp, tmp2, dfdu, dfdu_approx, dfdv, dfdv_approx);
        panic!();
    }
}
