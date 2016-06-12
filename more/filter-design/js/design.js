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

function findLowpassPolesAndZeroes(fp, fs, Ap, As, type) {
    var Wp = 2 * Math.PI * fp;
    var Ws = 2 * Math.PI * fs;
    var ep = Math.sqrt(Math.pow(10, Ap / 10) - 1);
    var es = Math.sqrt(Math.pow(10, As / 10) - 1);

    // Selectivity and discrimination parameters
    var k = Wp / Ws;
    var k1 = ep / es;

    var N;
    var K;
    var Kp;
    var K1;
    var K1p;

    // Determine order of filter to meet or exceed the requirements
    if (type === "butterworth") {
        N = Math.ceil(Math.log(1 / k1) / Math.log(1 / k));
    } else if (type === "cheby-1") {
        N = Math.ceil(Math.acosh(1 / k1) / Math.acosh(1 / k));
    } else if (type === "cheby-2") {
        N = Math.ceil(Math.acosh(1 / k1) / Math.acosh(1 / k));
        // Revompute k to satisfy degree equation
        k = 1 / Math.cosh(Math.acosh(1 / k1) / N);
    } else if (type === "elliptic") {
        K = elliptic_kc(k * k);
        Kp = elliptic_kc(1 - k * k);
        console.log("K = " + K + ", K' = " + Kp);

        K1 = elliptic_kc(k1 * k1);
        K1p = elliptic_kc(1 - k1 * k1);
        console.log("K1 = " + K1 + ", K1' = " + K1p);

        N = (K1p / K1) / (Kp / K);
        console.log("N = " + N + ", " + Math.ceil(N));
        N = Math.ceil(N);

        k = ellipticDeg(N, K1, K1p);
        K = elliptic_kc(k * k);
        console.log("new k = " + k + ", K = " + K);
    }

    var L = Math.floor(N / 2);
    var r = N - 2 * L;
    var u = new Array(L);
    var pa = new Array(L);
    var za;

    for (var m = 0; m < L; ++m) {
        u[m] = (2 * m + 1) / N;
    }

    // Determine poles and zeroes, if any
    if (type === "butterworth") {
        za = [];
        var factor = Ws / Math.pow(es, 1 / N);
        for (var m = 0; m < L; ++m) {
            var pr = factor * Math.cos(Math.PI / 2 * u[m]);
            var pi = factor * Math.sin(Math.PI / 2 * u[m]);
            pa[m] = makeComplex(-pi, pr);
            pa0 = -factor;
        }
    } else if (type === "cheby-1") {
        za = [];
        var v0 = Math.asinh(1 / ep) / (N * Math.PI / 2);
        for (var m = 0; m < L; ++m) {
            var pr = Math.cos(Math.PI / 2 * u[m]) * Math.cosh(Math.PI / 2 * v0);
            var pi = Math.sin(Math.PI / 2 * u[m]) * Math.sinh(Math.PI / 2 * v0);
            pa[m] = makeComplex(-Wp * pi, Wp * pr);
        }
        pa0 = -Wp * Math.sinh(v0 * Math.PI / 2);
    } else if (type === "cheby-2") {
        var factor = Wp / k;
        var v0 = Math.asinh(es) / (N * Math.PI / 2);
        var za = new Array(L);
        for (var m = 0; m < L; ++m) {
            za[m] = makeComplex(0, -factor / Math.cos(Math.PI / 2 * u[m]));
            pa[m] = cdiv(
                makeComplex(0, -factor),
                complex_cos(rcmul(Math.PI / 2, makeComplex(u[m], v0))));
        }
        pa0 = -factor / Math.sinh(v0 * Math.PI / 2);
    } else {
        // Elliptic
        var fs_new = fp / k;
        console.log("fs_new = " + fs_new);
        var L = Math.floor(N / 2);
        var r = N - 2 * L;

        var zeta_i = new Array(L);
        for (var n = 1; n <= L; ++n) {
            zeta_i[n - 1] = jacobi_cd((2 * n - 1) / N * K, k * k);
        }
        console.log("zeta_i = ");
        console.log(zeta_i);

        var za = zeta_i.map(function (z) {
            return makeComplex(0, Wp / (k * z));
        });

        console.log("za = ");
        console.log(za);

        var v0 = complex_inverse_jacobi_sn(makeComplex(0, 1 / ep), k1 * k1);
        v0 = v0.im / N / K1;
        console.log("v0 = ");
        console.log(v0);

        var pa0;
        pa0 = complex_jacobi_sn(makeComplex(0, v0 * K), makeComplex(k * k));
        pa0 = cmul(makeComplex(0, Wp) , pa0);
	pa0 = pa0.re;
        console.log("pa0 =");
        console.log(pa0);

        var pa = new Array(L);
        for (var n = 1; n <= L; ++n) {
            pa[n - 1] = complex_jacobi_cd(
                makeComplex((2 * n - 1) / N * K, -v0 * K),
                makeComplex(k * k));
            pa[n - 1] = cmul(makeComplex(0, Wp), pa[n - 1]);
        }
        console.log("pa = ");
        console.log(pa);
    }

    var H0 = 1;
    if (type === "cheby-1" || type === "elliptic") {
        H0 = Math.pow(Math.pow(10, -Ap / 20), 1 - r);
    }

    return {
        order: N,
        zeroes: za,
        poles: [pa0, pa],
        H0: H0
    }
}

function ellipticDeg(N, K1, K1p) {
    var q = Math.exp(-Math.PI*(K1p/K1)/N);
    // Compute the two sums
    //
    //   sum(q^(m*(m+1)), m, 0, 10) = 1 + sum(q^(m*(m+1)), m, 1, 10)
    //   sum(q^(m*m), m, 1, 10)
    var s1 = 1;
    var s2 = 0;
    for (var m = 1; m <= 10; ++m) {
        s1 += Math.pow(q, m*(m+1));
        s2 += Math.pow(q, m*m);
    }

    return 4*Math.sqrt(q)* Math.pow(s1 / (1 + 2*s2), 2);
}

function analogLowpassFilter(fp, fs, Ap, As, type) {
    var polesZeroes = findLowpassPolesAndZeroes(fp, fs, Ap, As, type);
    var N = polesZeroes.order;
    var za = polesZeroes.zeroes;
    var pa0 = polesZeroes.poles[0];
    var pa = polesZeroes.poles[1];
    var L = Math.floor(N / 2);
    var r = N - 2 * L;

    var A = [];
    var B = [];

    if (r === 1) {
        A.push([1, -1 / pa0]);
        B.push([1, 0]);
    }

    for (var m = 0; m < L; ++m) {
        var recip = crecip(pa[m]);
        A.push([1, -2 * recip.re, Math.pow(cabs(recip), 2)]);
    }

    if (type === "cheby-2" || type === "elliptic") {
        for (var m = 0; m < L; ++m) {
            var recip = crecip(za[m]);
            B.push([1, -2 * recip.re, Math.pow(cabs(recip), 2)]);
        }
    } else {
        for (var m = 0; m < L; ++m) {
            B.push([1, 0, 0]);
        }
    }

    return {
        order: N,
        H0: polesZeroes.H0,
        top: B,
        bot: A
    };
}

function texifyNumber(number, options) {
    // Convert the number to a string.  If it is in scientific form,
    // replace with the appropriate teX version.
    if (options && !options.showUnity && number === 1) {
        if (options.addSign) {
            return Math.sign(number) < 0 ? "-" : "+";
        }
        return "";
    }
    if (options && options.addSign) {
        var s = Math.abs(number).toString();
        s = s.replace(/(e)(.*)/, "\\times 10^{$2}");
        s = (number < 0) ? "- " + s : "+" + s;
    } else {
        var s = number.toString();
        s = s.replace(/(e)(.*)/, "\\times 10^{$2}");
    }

    return s;
}

function analogTermTeX(term) {
    var f = "\\frac{";
    if (term[0].length == 2) {
        // Linear term
        if (term[0][1] == 0) {
            f += "1";
        } else {
            f += "1 " + texifyNumber(term[0][1], {
                addSign: true
            }) + "s";
        }
        f += "}{";
        f += "1 " + texifyNumber(term[1][1], {
            addSign: true
        }) + "s";
        f += "}";
    } else {
        f += texifyNumber(term[0][0], {
            showUnity: true
        });
        if (term[0][1] != 0) {
            f += texifyNumber(term[0][1], {
                addSign: true
            }) + "s";
        }
        if (term[0][2] != 0) {
            f += texifyNumber(term[0][2], {
                addSign: true
            }) + "s^2";
        }
        f += "}{";
        f += "1 " + texifyNumber(term[1][1], {
            addSign: true
        }) + "s ";
        f += texifyNumber(term[1][2], {
            addSign: true
        }) + "s^2";
        f += "}";
    }

    return f;
}

function analogTeX(filter) {
    var f = "\\begin{align*}\n";
    f += "H_a(s) = \\, &";
    if (filter.H0 != 1) {
        f += texifyNumber(filter.H0) + "\\\\\n";
        f += "& \\times ";
    }
    f += analogTermTeX([filter.top[0], filter.bot[0]]) + "\\\\\n";
    for (var k = 1; k < filter.top.length; ++k) {
        var term = filter[k];
        f += "& \\times" + analogTermTeX([filter.top[k], filter.bot[k]]) + "\\\\\n";
    }
    f += "\n\\end{align*}\n";

    return f;
}

function analogResponse(filter, freq) {
    var top = filter.top;
    var bot = filter.bot;
    // To compute the response we need to calculate expression of the
    // form 1+b*s+c*s^2 where s = 2*pi*j*f.
    //
    //  1+b*s+c*s^2 = (1-4*pi*c*f^2) + i*2*pi*b*f
    var mag = new Float32Array(freq.length);
    var phase = new Float32Array(freq.length);
    mag.fill(filter.H0);
    phase.fill(0);

    var pi2 = 2 * Math.PI;
    var pi4 = pi2 * pi2;

    for (var k = 0; k < freq.length; ++k) {
        var f = freq[k];
        for (var m = 0; m < top.length; ++m) {
            if (top[m].length == 2) {
                mag[k] *= Math.hypot(1, top[m][1] * pi2 * f);
                mag[k] /= Math.hypot(1, bot[m][1] * pi2 * f);
            } else {
                // 1-4*pi*c*f^2 + i*2*pi*b*f.
                mag[k] *= Math.hypot(1 - top[m][2] * pi4 * f * f, pi2 * top[m][1] * f);
                mag[k] /= Math.hypot(1 - bot[m][2] * pi4 * f * f, pi2 * bot[m][1] * f);
                phase[k] += Math.atan2(pi2 * top[m][1] * f, 1 - top[m][2] * pi4 * f * f);
                phase[k] -= Math.atan2(pi2 * bot[m][1] * f, 1 - bot[m][2] * pi4 * f * f);
            }
        }
    }

    return {
        mag: mag,
        phase: phase
    };
}

function digitalLowpassFilter(fp, fs, Ap, As, Fs, type) {
    var wPass = 2 * Math.PI * fp / Fs;
    var wStop = 2 * Math.PI * fs / Fs;
    var omegaPass = Math.tan(wPass / 2);
    var omegaStop = Math.tan(wStop / 2);
    var polesZeroes = findLowpassPolesAndZeroes(omegaPass / (2 * Math.PI), omegaStop / (2 * Math.PI), Ap, As, type);
    var N = polesZeroes.order;
    var za = polesZeroes.zeroes;
    var pa0 = polesZeroes.poles[0];
    var pa = polesZeroes.poles[1];
    var H0 = polesZeroes.H0;

    var p0 = (1 + pa0) / (1 - pa0);
    var z = za.map(zz => cdiv(
        makeComplex(1 + zz.re, zz.im),
        makeComplex(1 - zz.re, -zz.im)));
    var p = pa.map(pp => cdiv(
        makeComplex(1 + pp.re, pp.im),
        makeComplex(1 - pp.re, -pp.im)));
    var G0 = (1 - p0) / 2;

    /*
        console.log("p0 = " + p0);
        console.log("z = ");
        console.log(z);
        console.log("pa =");
        console.log(pa);
        console.log("p = ");
        console.log(p);
    */

    var A = [];
    var B = [];

    if ((N % 2) == 1) {
        A.push([1, -p0]);
        B.push([G0, [1, 1]]);
    }

    if (type === "butterworth" || type === "cheby-1") {
        var G = p.map(pp => {
            return makeComplex((1 - pp.re) / 2, -pp.im / 2)
        });
        console.log("G =");
        console.log(G);
        for (var m = 0; m < p.length; ++m) {
            var g = Math.pow(cabs(G[m]), 2);
            A.push([1, -2 * p[m].re, Math.pow(cabs(p[m]), 2)]);
            B.push([g, [1, 2, 1]]);
        }
    } else if (type === "cheby-2" || type === "elliptic") {
        for (var m = 0; m < p.length; ++m) {
            var G = cdiv(
                makeComplex(1 - p[m].re, -p[m].im),
                makeComplex(1 - z[m].re, -z[m].im));
            var g = Math.pow(cabs(G), 2);
            A.push([1, -2 * p[m].re, Math.pow(cabs(p[m]), 2)]);
            // For cheby-2 filters, I think we have notch filters, so
            // the coefficient of the z^(-2) term should be exactly 0.
            // Make it so.
            var c = (type === "cheby-2") ? 1 : Math.pow(cabs(z[m]), 2);
            B.push([g, [1, -2 * z[m].re, c]]);
        }
    }

    return {
        order: N,
        H0: H0,
        top: B,
        bot: A
    };
}

function digitalTermTeX(term) {
    var f = "\\frac{";
    if (term[0][1].length == 2) {
        // Linear term
        if (term[0][1][1] == 0) {
            f += "1";
        } else {
            f += "1 " + texifyNumber(term[0][1][1], {
                addSign: true
            }) + "z^{-1}";
        }
        f += "}{";
        f += "1 " + texifyNumber(term[1][1], {
            addSign: true
        }) + "z^{-1}";
        f += "}";
    } else {
        f += texifyNumber(term[0][0], {
            showUnity: true
        }) + "(1";
        if (term[0][1][0] != 0) {
            f += texifyNumber(term[0][1][1], {
                addSign: true
            }) + "z^{-1}";
        }
        if (term[0][1][2] != 0) {
            f += texifyNumber(term[0][1][2], {
                addSign: true
            }) + "z^{-2}";
        }
        f += ")}{";
        f += "1 " + texifyNumber(term[1][1], {
            addSign: true
        }) + "z^{-1} ";
        f += texifyNumber(term[1][2], {
            addSign: true
        }) + "z^{-2}";
        f += "}";
    }

    return f;
}

function digitalTeX(filter) {
    var f = "\\begin{align*}\n";
    f += "H(z) = \\, &";

    if (filter.H0 != 1) {
        f += texifyNumber(filter.H0) + "\\\\\n";
        f += "&\\times";
    }
    f += digitalTermTeX([filter.top[0], filter.bot[0]]) + "\\\\\n";
    for (var k = 1; k < filter.top.length; ++k) {
        var term = filter[k];
        f += "& \\times" + digitalTermTeX([filter.top[k], filter.bot[k]]) + "\\\\\n";
    }
    f += "\n\\end{align*}\n";

    return f;
}

function digitalResponse(filter, freq, Fs) {
    // To compute the response of the digital filter we need to evaluate terms like
    // 1 + b*z^(-1)+c*z^(-2) for z = exp(j*w).
    //
    //  1+b*z^(-1)+c*z^(-2) = (1+b*cos(w)+c*cos(2*w)) - j*(b*sin(w)+c*sin(2*w)
    var top = filter.top;
    var bot = filter.bot;
    var mag = new Float32Array(freq.length);
    var phase = new Float32Array(freq.length);
    mag.fill(filter.H0);

    for (var k = 0; k < freq.length; ++k) {
        var w = 2 * Math.PI * freq[k] / Fs;
        for (var m = 0; m < top.length; ++m) {
            var g = top[m][0];
            var B = top[m][1];
            var A = bot[m];
            var sumt = makeComplex(1);
            var sumb = makeComplex(1);
            for (var n = 1; n < B.length; ++n) {
                var cis = makeComplex(Math.cos(n * w), -Math.sin(n * w));
                sumt = cadd(sumt, rcmul(B[n], cis));
                sumb = cadd(sumb, rcmul(A[n], cis));
            }
            mag[k] *= g * cabs(sumt) / cabs(sumb);
            phase[k] += Math.atan2(sumt.im, sumt.re) - Math.atan2(sumb.im, sumb.re);
        }
    }

    return {
        mag: mag,
        phase: phase
    };
}

function webAudioFilterDesc(top, bot, Fs, type) {
    var order = bot.length - 1;
    var gain = top[0];
    var zterm = top[1];
    if (order == 1) {
        return {
            filterType: "iir",
            top: [gain, gain],
            bot: bot,
            filterGain: 1
        };
    }

    if (type === "butterworth" || type === "cheby-1") {
        var b = bot[1];
        var c = bot[2];
        var alpha = (1 - c) / (1 + c);
        var w0 = Math.acos(-b * (1 + alpha) / 2);
        var a0 = 1 + alpha;
        var b0 = (1 - Math.cos(w0)) / 2;
        return {
            filterType: "biquad",
            biquadType: "lowpass",
            gain: gain * a0 / b0,
            f0: w0 * Fs / (2 * Math.PI),
            Q: 20 * Math.log10(Math.sin(w0) / 2 / alpha),
            top: zterm,
            bot: bot,
            filterGain: gain
        };
    }
    if (type === "cheby-2" || type === "elliptic") {
        // Superficially, sections for a Chebyshev-2 and elliptic
        // filters kind of look like biquad notch filters, but they're
        // not.  The numerator coefficient of z^(-1) is not consistent
        // with the denominator coefficient of z^(-1}.  And sometimes
        // the numerator coefficient is positive, whereas the biquad
        // notch filter has a negative coefficient.
        //
        // Thus, use an IIR filter.
        return {
            filterType: "iir",
            top: zterm.map(x => x * gain),
            bot: bot,
            filterGain: 1
        }
    }
    throw "Unknown filter type: " + type;
}

function webAudioFilter(filter, Fs, type) {
    var result = [];

    for (var m = 0; m < filter.top.length; ++m) {
        result.push(webAudioFilterDesc(filter.top[m], filter.bot[m], Fs, type));
    }

    var totalGain = filter.H0;
    for (var m = 0; m < result.length; ++m) {
        if (result[m].filterType === "biquad") {
            totalGain *= result[m].gain;
        }
    }

    return {
        totalGain: totalGain,
        desc: result
    }
}
