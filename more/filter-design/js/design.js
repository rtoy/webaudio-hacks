// Based on Sophocles J. Orfandidis, "Lecture Notes on Elliptic Filter
// Design", Rutgers University.
// www.ece.rutgers.edu/~orfanidi/ece521/notes.pdf
//
//
// The magnitude response for an analog lowpass Butterworth,
// Chebyshev, or elliptic filter is
//
//   |H(omega)|^2 = 1/(1+eps_p^2*F(N,w)^2), w = omega/omega_p
//
// where N is the filter order, F(N,w) is a function of the normalized
// frequency w given by:
//
// Butterworth:  w^N
// Chebyshev-1:  C(N,w)
// Chebyshev-2:  1/[k1*C(N, 1/k/w)]
// Elliptic:     cd(N*u*K1, k1), w = cd(u*K, k)
//
// where C(N, x) is the order-N Chebyshev polynomial, hat is C(N, x) =
// cos(N*acos(x)). and cd(x,k) is the Jacobian elliptic function cd
// with modulus k and real quarter-period K.
//
// The typical filter design calls for omega_p, the passband edge,
// omega_s, the stopband edge, with the corresponding gains G_p and
// G_s.  These can be determined from the attenuations A_p and
// A_s. A_p is the maximum allowed deviation of the passband (the
// passband should lie between 0 and -A_p dB) and A_s is the minimum
// attentuation of the stop band (the stopband gain should always be
// less than -A_s dB).
//
// The ripple parameters eps_p and eps_s are
//
//   G_p = 1/sqrt(1+eps_p^2) = 10^(-A_p/20)
//   G_s = 1/sqrt(1+eps_s^2) = 10^(-A_s/20)
//
// or equivalently
//
//   A_p = -20*log10(G_p) = 10*log10(1+eps_p^2)
//   A_s = -20*log10(G_s) = 10*log10(1+eps_s^2)
//   eps_p = sqrt(1/G_p^2-1) = sqrt(10^(A_p/10)-1)
//   eps_s = sqrt(1/G_s^2-1) = sqrt(10^(A_s/10)-1)
//
// Define the following design parameters:
//
//   k = omega_p / omega_s
//   k1 = eps_p / eps_s
//
// k is known as the selectivity and k1 is the discrimination. k <= 1
// and k1 < 1.
//
// The requirement that the passband and stopband specifications are
// met at the corners omega_p (w = 1) and omega_s (w = 1/k) gives rise
// to the conditions
//
//   |H(omega_p)|^2 = 1/(1+eps_p^2*F(N,1)^2) = 1/(1+eps_p^2) => F(N,1) = 1
//   |H(omega_s)|^2 = 1/(1+eps_s^2*F(N,1/k)^2) = 1/(1+eps_s^2) => F(N,1/k) = 1/k1
//
// Thus, F(N,w) is normalized such that F(N,1)=1 and must satisfy the
// degree equation that relates N, k, and k1:
//
//   F(N, 1/k) = 1/k1.
//
// For a Butterworth filter:
//
//   k^(-N) = 1/k1 => N = log(1/k1)/log(1/k)
//
// For both Chebyshev cases:
//
//   C(N,1/k) = 1/k1 => N = acosh(1/k1)/acosh(1/k)
//
// In practice, one specifies omega_p, omega_s, eps_p and eps_s which
// fix the values of k and k1.  Thus, we solve for N to determine the
// filter orde, which must be rounded up to the next integer.  This
// necessitates recomputing either k in terms of N and k1 or k1 in
// terms of N and k.
//
// For the elliptic case, we want a filter that is equiripple in both
// the passband and stopband.  This is achieved if
//
//   F(N, w) = 1/(k1*F(N, 1/(k1*w)))
//
// The degree equation is then
//
//   N*K'/K = K1'/K1
//
// where K, K1 are the complete elliptic integrals corresponding to
// the moduli k, k1, and K', K1' are the complete elliptic integrals
// corresponding to the complementary moduli k' = sqrt(1-k^2) and k1'
// = sqrt(1-k1^2).
//
// The degree equation can be written as
//
//   k1 = k^N*product(sn(u_i*K, k)^4), i = 1, L)
//
// Where N = 2*L + r, u_i = (2*i-1)/N
//
//----------------------------------------------------------------------
// Analog Filter Design
//
// The transfer function is constructed from the zeros and poles:
//
//   H_a(s) = H0*[1/(s-p[0])]^r
//     * product((1-s/z[i])*(1-s/conj(z[i]))/(1-s/p[i])/(1-s/conj(p[i])), i, 1, L)
//
// where N = 2*L + r, r = 0, 1, and z[i] and p[i] are the poles and
// zeroes of the filter.
//
// H0 is the DC gain and H0 = 1 for Butterworth and Chebyshev-2
// filters and G_p^(1-r) for Chebyshev-1 and elliptic filters.
//
// Equivalently:
//
//   H(s) = H0*[1/(1+A[0]*s)]^5
//    * product((1+B1[i]*s+B2[i]*s^2)/(1+A1[i]*s+A2[i]*s^2),i,1,L)
//
// where
//   [1, B1[i], B2[i]] = [1, -2*Re(1/z[i]), 1/|z[i]|^2]
//   [1, A1[i], A2[i]] = [1, -2*Re(1/p[i]), 1/|p[i]|^2]
//   [1, A0[0]] = [1, -1/p[0]]
//
// The zeroes z[i] will arise from the poles of F(N,w) and the poles
// p[i] will arise from the zeroes of 1+eps_p^2*F(N,w)^2 = 0.
//
// The poles of F(N,w) occur at w[i] = 1/(k*zeta[i])m, zeta[i] =
// cd(u_i*K, k).
//
// The poles p[i] are found by solving
//
//   F(N, w) = +/- j/eps_p.
//
// The left-hand poles are then
//
//   p[i] = omega_p*j*w*cd((u_i - j*v_0)*K,k)
//
// where v0 is the real-valued solution of
//
//   sn(j*v0*N*K1, k1) = j/eps_p
//
// or
//
//   v0 = -j/(N*K1)*inverse_sn(j/eps_p, k1)
//
//
//----------------------------------------------------------------------
// Jacobian elliptic fucntions.
//
function jacobi_sn(u, m) {
    if (m == 0) {
        // sn(u, 0) = sin(u)
        return Math.sin(u);
    }
    if (m == 1) {
        // sn(u, 1) = tanh(u)
        return Math.tanh(u);
    }
    return  elliptic_sn_descending(u, m);
}

function jacobi_snk(u, k) {
    return jacobi_sn(u, k * k);
}

function elliptic_sn_descending(u, m) {
    if (m == 1) {
        // sn(u, 1) = tanh(u)
        return Math.tanh(u);
    }
    if (Math.abs(m) <= epsilon*Math.abs(u)) {
        return Math.sinh(u);
    }

    var [v, mu, root_mu] = descending_transform(u, m);
    var new_sn = elliptic_sn_descending(v, mu);
    return (1 + root_mu) * new_sn / (1 + root_mu * new_sn * new_sn);
}

// Descending Landen transform
function descending_transform(u, m) {
    var root_m1 = Math.sqrt(1 - m);
    var root_mu = (1 - root_m1) / (1 + root_m1);
    var mu = root_mu * root_mu;
    var v = u / (1 + root_mu);

    return {v: v, mu: mu, root_mu: root_mu};
}
        
function jacobi_cn(u, m) {
    if (m == 0) {
        // cn(u, 0) = cos(u)
        return Math.cos(u);
    }
    if (m == 1) {
        // cn(u, 1) = sech(u)
        return 1 / Math.cosh(u);
    }
    var [v, mu, root_mu] = ascending_transform(u, m);
    var d = jacobi_dn(v, mu);
    return (1 + root_mu1) / mu * ((d * d - root_mu1) / d);
}

function jacobi_cnk(u, k) {
    return jacobi_cn(u, k*k);
}

function ascending_transform(u, m) {
    var root_m = Math.sqrt(m);
    var mu = 4 * root_m / Math.pow(1 + root_m, 2);
    var root_mu1 = (1 - root_m) / (1 + root_m);
    var v = u / (1 + root_mu1);

    return {v: v, mu: mu, root_mu1: root_mu1};
}

function jacobi_dn(u, m) {
    if (m == 0) {
        // dn(u, 0) = 1
        return 1;
    }
    if (m == 1) {
        // dn(u, 1) = sech(u)
        return 1 / Math.cosh(u);
    }
    var root_1m = Math.sqrt(1 - m);
    var root = (1 - root_1m) / (1 + root_1m);
    var z = u / (1 + root);
    var s = elliptic_sn_ascending(z, root * root);
    var p = root * s * s;
    return (1 - p) / (1 + p);
}

function jacobi_dnk(u, k) {
    return jacobi_dn(u, k*k);
}

funtion jacobi_cdk(u, k) {
    return jacobi_cnk(u, k) / jacobi_dnk(u, k);
}

function elliptic_kc(m) {
    if (< m 0) {
        m = -m;
        var m1 = 1 + m;
        var root = Math.sqrt(m1);
        var mdiv = m / m1;
        return elliptic_kc(mdiv) / root - (elliptic_f 0 mdiv) / root;
    }
    if (m == 0) {
        return Math.PI / 2;
    }
    if (m == 1) {
        throw ;
    }
    return carlson_rf(0, 1 - m, 1);
}

function elliptic_kck(k) {
    return elliptic_kc(k*k);
}

function carlson_rf(x, y, z) {
    var xn = x;
    var yn = y;
    var zn = z;
    var a = (xn + yn + zn) / 3;
    var eps = Math.max(Math.abs(a - xn), Math.abs(a - yn), Math.abs(a - zn)) / errtol(x, y, z);
    var an = a;
    var power4 = 1;
    var n = 0;

    while (power4 * eps > Math.abs(an)) {
        var xnroot = Math.sqrt(xn);
        var ynroot = Math.sqrt(yn);
        var znroot = Math.sqrt(zn);
        var lam = xnroot * xnroot + ynroot * ynroot + znroot * znroot;
        power4 *= 1/4;
        xn = (xn + lam) / 4;
        yn = (yn + lam) / 4;
        zn = (zn + lam) / 4;
        an = (an + lam) / 4;
        ++n;
    }
    var xndev = (a - x) * power4 / an;
    var yndev = (a - y) * power4 / an;
    var zndev = -(xndev + yndev);
    var ee2 = xndev*yndev - 6*zndev*zndev;
    var ee3 = xndev * yndev * zndev;
    var s = 1 - ee2/10 + ee3/14 + ee2*ee2/24 - 3/44*ee2*ee3;
    return s / Math.sqrt(an);
}

function inverse_jacobi_sn(u, m) {
    return u * carlson_rf(1 - u*u, 1 - m*u*u, 1);
}

function inverse_jacobi_snk(u, k) {
    return inverse_jacobi_sn(u, k*k);
}

function inverse_jacobi_dn(w, m) {
    if (w == 1) {
        return 0;
    }
    if (m == 1) {
        // w = dn(u,1) = sech(u) = 1 / cosh(u)
        // cosh(u) = 1/w
        // u = acosh(1/w)
        return Math.acosh(1/w;
    }
    return inverse_jacobi_cn(w, 1/m) / Math.sqrt(m);
        
}

function inverse_jacobi_cd(w, m) {
    return inverse_jacobi_sn(Math.sqrt(1 - w*w) / Math.sqrt(1 - m*w*w), m);
}

function inverse_jacobi_cdk(w, k) {
    return inverse_jacobi_cd(w, k*k);
}
