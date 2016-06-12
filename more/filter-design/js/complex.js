// Very rudimentary complex math package that does just enough for what I want.
//
// A complex number is represented as a dictionary: {re: r, im: i}.

// Make a complex number.  |y| is optional.
function makeComplex(x, y) {
    return {
        re: x,
        im: y || 0
    };
}

// Compare two complex numbers x and y for equality.  Returns true if
// the real and imaginary components are equal.
function cequal(x, y) {
    return (x.re == y.re) && (x.im == y.im);
}

// Compare a complex number \x| against a real number |r|, returning
// true if the real part of |x| is equal to r and the imaginary part
// is 0.
function cequalr(x, r) {
    return (x.re == r) && (x.im == 0);
}

// Complex addition: x + y
function cadd(x, y) {
    return makeComplex(x.re + y.re, x.im + y.im);
}

// Complex addition of real and complex: r + z
function rcadd(r, z) {
    return makeComplex(r + z.re, z.im);
}

// Complex subtraction: x - y
function csub(x, y) {
    return makeComplex(x.re - y.re, x.im - y.im);
}

// Complex subtraction of real and complex: r - z
function rcsub(r, z) {
    return makeComplex(r - z.re, -z.im);
}

// Complex multiplication: x*y
function cmul(x, y) {
    return makeComplex(
        x.re * y.re - x.im * y.im,
        x.im * y.re + x.re * y.im
    );
}

// Reciprocal of complex number z
function crecip(z) {
    var d = z.re * z.re + z.im * z.im;
    return makeComplex(z.re / d, -z.im / d);
}

// Absolute value of a complex number z
function cabs(z) {
    return Math.hypot(z.re, z.im);
}

// Multiply a real times a complex: r*x
function rcmul(a, x) {
    return makeComplex(a * x.re, a * x.im);
}
    
// Complex division: w/z
function cdiv(w, z) {
    var den = Math.pow(Math.hypot(z.re, z.im), 2);
    return makeComplex(
        (w.re * z.re + w.im * z.im) / den,
        (w.im * z.re - w.re * z.im) / den
    );
}

// Divide complex by real
function cdivr(w, r) {
    return makeComplex(w.re / r, w.im / r);
}

// Divide real by complex
function rcdiv(r, z) {
    var den = Math.pow(Math.hypot(z.re, z.im), 2);
    return makeComplex(z.re / den, -z.im / den);
}

// Principal square root of a complex number z.
function csqrt(z) {
    var m = cabs(z);
    var r = Math.sqrt((m + z.re) / 2);
    var i = Math.sqrt((m - z.re) / 2);
    return makeComplex(r, Math.sign(z.im) * i);
}

// Sin of a complex value.
function complex_sin(z) {
    // sin(x+%i*y) = sin(x)*cosh(y) + i*cos(x)*sinh(y)
    // For our purposes, this is good enough.
    return makeComplex(
        Math.sin(z.re) * Math.cosh(z.im),
        Math.cos(z.re) * Math.sinh(z.im)
    );
}

// Cos of a complex value
function complex_cos(u) {
    // cos(x+%i*y) = cos(x)*cosh(y) - i*sin(x)*sinh(y)
    // For our purposes, this is good enough.
    return makeComplex(
        Math.cos(z.re) * Math.cosh(z.im),
	-Math.sin(z.re) * Math.sinh(z.im)
    );
}

// Cosh of a complex value
function complex_cosh(z) {
    // cosh(x+i*y) = cosh(x)*cos(y) + i*sinh(x)*sin(y)
    return makeComplex(
        Math.cosh(z.re) * Math.cos(z.im),
        Math.snih(z.re) * Math.cos(z.im)
    );
}

// Tanh of a complex value
function complex_tanh(z) {
    // tanh(x+i*y) = (tanh(x)+i*%tan(y))/(1+i*tanh(x)*tan(y))
    // For our purposes, this is good enough.
    return cdiv(
	makeComplex(Math.tanh(z.re), Math.tan(z.im)),
        makeComplex(1, Math.tanh(z.re) * Math.tan(z.im))
    );
}
