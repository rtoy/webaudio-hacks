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
    var Wp = 2*Math.PI*fp;
    var Ws = 2*Math.PI*fs;
    var ep = Math.sqrt(Math.pow(10, Ap/10)-1);
    var es = Math.sqrt(Math.pow(10, As/10)-1);

    // Selectivity and discrimination parameters
    var k = Wp/Ws;
    var k1 = ep/es;

    var N;

    // Determine order of filter to meet or exceed the requirements
    if (type === "butterworth") {
	N = Math.ceil(Math.log(1/k1)/Math.log(1/k));
    } else if (type === "cheby-1") {
	N = Math.ceil(Math.acosh(1/k1)/Math.acosh(1/k));
    } else if (type === "cheby-2") {
	N = Math.ceil(Math.acosh(1/k1)/Math.acosh(1/k));
	// Revompute k to satisfy degree equation
	k = 1/Math.cosh(Math.acosh(1/k1) / N);
    } else if (type === "elliptic") {
	throw "Unknown filter type: " + type;
    }

    var L = Math.floor(N/2);
    var r = N - 2*L;
    var u = new Array(L);
    var pa = new Array(L);
    var za;

    for (var m = 0; m < L; ++m) {
	u[m] = (2*m+1)/N;
    }

    // Determine poles and zeroes, if any
    if (type === "butterworth") {
	za = [];
	var factor = Ws / Math.pow(es, 1/N);
	for (var m = 0; m < L; ++m) {
	    var pr = factor*Math.cos(Math.PI/2*u[m]);
	    var pi = factor*Math.sin(Math.PI/2*u[m]);
	    pa[m] = {re: -pi, im: pr};
	    pa0 = -factor;
	}
    } else if (type === "cheby-1") {
	za = [];
	var v0 = Math.asinh(1/ep)/(N*Math.PI/2);
	for (var m = 0; m < L; ++m) {
	    var pr = Math.cos(Math.PI/2*u[m])*Math.cosh(Math.PI/2*v0);
	    var pi = Math.sin(Math.PI/2*u[m])*Math.sinh(Math.PI/2*v0);
	    pa[m] = {re: -Wp*pi, im: Wp*pr};
	}
	pa0 = -Wp*Math.sinh(v0*Math.PI/2);
    } else if (type === "cheby-2") {
	var factor = Wp / k;
	var v0 = Math.asinh(es) / (N*Math.PI/2);
	var za = new Array(L);
	for (var m = 0; m < L; ++m) {
	    za[m] = {re: 0, im: -factor/Math.cos(Math.PI/2*u[m])};
	    var su = Math.sin(Math.PI/2*u[m]);
	    var cu = Math.cos(Math.PI/2*u[m]);
	    var shv = Math.sinh(Math.PI/2*v0);
	    var chv = Math.cosh(Math.PI/2*v0);
	    var d = su*su*shv*shv + cu*cu*chv*chv;
	    var pr = -factor*su*shv/d;
	    var pi = -factor*cu*chv/d;
	    pa[m] = {re: pr, im: pi};
	}
	pa0 = -factor / Math.sinh(v0*Math.PI/2);
    } else {
	// Elliptic
    }

    var H0 = 1;
    if (type === "cheby-1" || type === "elliptic") {
	H0 = Math.pow(Math.pow(10, -Ap/20), 1-r);
    }

    return {order: N, zeroes: za, poles: [pa0, pa], H0: H0}
}

// Reciprocal of complex number z
function crecip(z) {
    var d = z.re*z.re + z.im*z.im;
    return {re: z.re/d, im: -z.im/d};
}

// Absolute value of a complex number z
function cabs(z) {
    return Math.hypot(z.re, z.im);
}

// Complex division: w/z
function cdiv(w, z) {
    var den = Math.pow(Math.hypot(z.re, z.im), 2);
    return {re: (w.re*z.re + w.im*z.im)/den,
	    im: (w.im*z.re - w.re*z.im)/den};
}

function analogLowpassFilter(fp, fs, Ap, As, type) {
    var polesZeroes = findLowpassPolesAndZeroes(fp, fs, Ap, As, type);
    var N = polesZeroes.order;
    var za = polesZeroes.zeroes;
    var pa0 = polesZeroes.poles[0];
    var pa = polesZeroes.poles[1];
    var L = Math.floor(N/2);
    var r = N - 2*L;

    var A = [];
    var B = [];

    if (r === 1) {
	A.push([1, -1/pa0, 0]);
	B.push([1,0,0]);
    }

    for (var m = 0; m < L; ++m) {
	var recip = crecip(pa[m]);
	A.push([1, -2*recip.re, Math.pow(cabs(recip),2)]);
    }

    if (type === "cheby-2" || type === "elliptic") {
	for (var m = 0; m < L; ++m) {
	    var recip = crecip(za[m]);
	    B.push([1, -2*recip.re, Math.pow(cabs(recip),2)]);
	}
    } else {
	for (var m = 0; m < L; ++m) {
	    B.push([1,0,0]);
	}
    }

    return {order: N, top: B, bot: A};
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
	s = (number < 0) ? "- " + s: "+" + s;
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
	    f += "1 " + texifyNumber(term[0][1], {addSign: true}) + "s";
	}
	f += "}{";
	f += "1 " + texifyNumber(term[1][1], {addSign: true}) + "s";
	f += "}";
    } else {
	f += texifyNumber(term[0][0], {showUnity: true});
	if (term[0][1] != 0) {
	    f += texifyNumber(term[0][1], {addSign: true}) + "s";
	}
	if (term[0][2] != 0) {
	    f += texifyNumber(term[0][2], {addSign: true}) + "s^2";
	}
	f += "}{";
	f += "1 " + texifyNumber(term[1][1], {addSign: true}) + "s ";
	f += texifyNumber(term[1][2], {addSign: true}) + "s^2";
	f += "}";
    }

    return f;
}

function analogTeX(filter) {
    var f = "\\begin{align*}\n";
    f += "H_a(s) = ";
    f += "& " + analogTermTeX([filter.top[0], filter.bot[0]]) + "\\\\\n";
    for (var k = 1; k < filter.top.length; ++k) {
	var term = filter[k];
	f += "& \\times" + analogTermTeX([filter.top[k], filter.bot[k]]) + "\\\\\n";
    }
    f += "\n\\end{align*}\n";

    return f;
}

function digitalLowpassFilter(fp, fs, Ap, As, Fs, type) {
    var wPass = 2*Math.PI*fp/Fs;
    var wStop = 2*Math.PI*fs/Fs;
    var omegaPass = Math.tan(wPass/2);
    var omegaStop = Math.tan(wStop/2);
    var polesZeroes = findLowpassPolesAndZeroes(omegaPass/(2*Math.PI), omegaStop/(2*Math.PI), Ap, As, type);
    var N = polesZeroes.order;
    var za = polesZeroes.zeroes;
    var pa0 = polesZeroes.poles[0];
    var pa = polesZeroes.poles[1];
    var H0 = polesZeroes.H0;
    
    var p0 = (1+pa0)/(1-pa0);
    var z = za.map(zz => cdiv({re: 1+zz.re, im: zz.im}, {re: 1-zz.re, im: -zz.im}));
    var p = pa.map(pp => cdiv({re: 1+pp.re, im: pp.im}, {re: 1-pp.re, im: -pp.im}));
    var G0 = (1 - p0)/2;

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
	var G = p.map(pp => { return {re: (1-pp.re)/2, im: -pp.im/2}; });
	console.log("G =");
	console.log(G);
	for (var m = 0; m < p.length; ++m) {
	    var g = Math.pow(cabs(G[m]), 2);
	    A.push([1, -2*p[m].re, Math.pow(cabs(p[m]),2)]);
	    B.push([g, [1, 2, 1]]);
	}
    } else if (type === "cheby-2" || type === "elliptic") {
	var G = p.map(pp => { return {re: (1-pp.re)/2, im: -pp.im/2}; });
	for (var m = 0; m < p.length; ++m) {
	    var g = Math.pow(cabs(G[m]), 2);
	    A.push([1, -2*p[m].re, Math.pow(cabs(p[m]),2)]);
	    B.push([g, [1, -2*z[m].re, Math.pow(cabs(z[m]),2)]]);
	}
    }

    return {order: N, H0: H0, top: B, bot: A};
}

function digitalTermTeX(term) {
    var f = "\\frac{";
    if (term[0][1].length == 2) {
	// Linear term
	if (term[0][1][1] == 0) {
	    f += "1";
	} else {
	    f += "1 " + texifyNumber(term[0][1][1], {addSign: true}) + "z^{-1}";
	}
	f += "}{";
	f += "1 " + texifyNumber(term[1][1], {addSign: true}) + "z^{-1}";
	f += "}";
    } else {
	f += texifyNumber(term[0][0], {showUnity: true}) + "(1";
	if (term[0][1][0] != 0) {
	    f += texifyNumber(term[0][1][1], {addSign: true}) + "z^{-1}";
	}
	if (term[0][1][2] != 0) {
	    f += texifyNumber(term[0][1][2], {addSign: true}) + "z^{-2}";
	}
	f += ")}{";
	f += "1 " + texifyNumber(term[1][1], {addSign: true}) + "z^{-1} ";
	f += texifyNumber(term[1][2], {addSign: true}) + "z^{-2}";
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
    f +=  digitalTermTeX([filter.top[0], filter.bot[0]]) + "\\\\\n";
    for (var k = 1; k < filter.top.length; ++k) {
	var term = filter[k];
	f += "& \\times" + digitalTermTeX([filter.top[k], filter.bot[k]]) + "\\\\\n";
    }
    f += "\n\\end{align*}\n";

    return f;
}

function webAudioFilterDesc(top, bot, Fs, type) {
    var order = bot.length - 1;
    var gain = top[0];
    var zterm = top[1];
    if (order == 1) {
	return {filterType: "iir",
		top: [gain, gain],
		bot: [1, zterm[0]]};
    }
    
    if (type === "butterworth" || type === "cheby-1") {
	var b = bot[1];
	var c = bot[2];
	var alpha = (1-c)/(1+c);
	var w0 = Math.acos(b*(1+alpha)/2);
	var a0 = 1+alpha;
	var b0 = (1-Math.cos(w0))/2;
	return {filterType: "biquad",
		biquadType: "lowpass",
		gain: gain*a0/b0,
		f0: w0*Fs/(2*Math.PI),
		Q: 20*Math.log10(Math.sin(w0)/2/alpha)};
    }
    if (type === "cheby-2" || type === "elliptic") {
	var b = bot[1];
	var c = bot[2];
	var alpha = (1-c)/(1+c);
	var w0 = Math.acos(b*(1+alpha)/2);
	var a0 = 1+alpha;
	return {filterType: "biquad",
		biquadType: "notch",
		gain: gain*a0,
		f0: w0*Fs/(2*Math.PI),
		Q: Math.sin(w0)/2/alpha};
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

    return {totalGain: totalGain,
	    desc: result}
}
