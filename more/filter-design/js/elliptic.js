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
    return elliptic_sn_descending(u, m);
}

function complex_jacobi_sn(u, m) {
    if (m == 0) {
        return complex_sin(u);
    }
    if (m == 1) {
        return complex_tanh(u);
    }
    return complex_elliptic_sn_descending(u, m);
}

function jacobi_snk(u, k) {
    return jacobi_sn(u, k * k);
}

function elliptic_sn_descending(u, m) {
    if (m == 1) {
        // sn(u, 1) = tanh(u)
        return Math.tanh(u);
    }
    if (Math.abs(m) <= Number.EPSILON * Math.abs(u)) {
        return Math.sin(u);
    }

    var {v, mu, root_mu} = descending_transform(u, m);
    var new_sn = elliptic_sn_descending(v, mu);
    return (1 + root_mu) * new_sn / (1 + root_mu * new_sn * new_sn);
}

function complex_elliptic_sn_descending(u, m) {
    if (m == 1) {
        // sn(u, 1) = tanh(u)
        return complex_tanh(u);
    }
    if (Math.abs(m) <= Number.EPSILON * Math.abs(u)) {
        return complex_sin(u);
    }

    var {v, mu, root_mu} = complex_descending_transform(u, m);
    var new_sn = complex_elliptic_sn_descending(v, mu);
    return cdiv(cmul(rcadd(1, root_mu), new_sn),
                rcadd(1, cmul(root_mu, cmul(new_sn, new_sn))));
}

// Descending Landen transform
function descending_transform(u, m) {
    var root_m1 = Math.sqrt(1 - m);
    var root_mu = (1 - root_m1) / (1 + root_m1);
    var mu = root_mu * root_mu;
    var v = u / (1 + root_mu);

    return {
        v: v,
        mu: mu,
        root_mu: root_mu
    };
}

function complex_descending_transform(u, m) {
    var root_m1 = csqrt(rcsub(1, m));
    var root_mu = cdiv(rcsub(1, root_m1), rcadd(1, root_m1));
    var mu = cmul(root_mu, root_mu);
    var v = cdiv(u, rcadd(1, root_mu));

    return {
        v: v,
        mu: mu,
        root_mu: root_mu
    };
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
    var {v, mu, root_mu1} = ascending_transform(u, m);
    var d = jacobi_dn(v, mu);
    return (1 + root_mu1) / mu * ((d * d - root_mu1) / d);
}

function complex_jacobi_cn(u, m) {
    if (m == 0) {
        return complex_cos(u);
    }
    if (m == 1) {
        return rcdiv(1, complex_cosh(u));
    }
    var {v, mu, root_mu1} = complex_ascending_transform(u, m);
    var d = complex_jacobi_dn(v, mu);
    return cmul(cdiv(rcadd(1, root_mu1), mu),
                cdiv(csub(cmul(d, d), root_mu1), d));
}

function jacobi_cnk(u, k) {
    return jacobi_cn(u, k * k);
}

function ascending_transform(u, m) {
    var root_m = Math.sqrt(m);
    var mu = 4 * root_m / Math.pow(1 + root_m, 2);
    var root_mu1 = (1 - root_m) / (1 + root_m);
    var v = u / (1 + root_mu1);

    return {
        v: v,
        mu: mu,
        root_mu1: root_mu1
    };
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
    var s = elliptic_sn_descending(z, root * root);
    var p = root * s * s;
    return (1 - p) / (1 + p);
}

function jacobi_dnk(u, k) {
    return jacobi_dn(u, k*k);
}

function jacobi_cdk(u, k) {
    return jacobi_cnk(u, k) / jacobi_dnk(u, k);
}

function elliptic_kc(m) {
    if (m < 0) {
        m = -m;
        var m1 = 1 + m;
        var root = Math.sqrt(m1);
        var mdiv = m / m1;
        return elliptic_kc(mdiv) / root - elliptic_f(0, mdiv) / root;
    }
    if (m == 0) {
        return Math.PI / 2;
    }
    if (m == 1) {
        return NaN;
    }
    return carlson_rf(0, 1 - m, 1);
}

function elliptic_kck(k) {
    return elliptic_kc(k*k);
}

function errtol(x, y, z) {
    return Math.sqrt(Number.EPSILON);
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
        var lam = xnroot * ynroot + xnroot * znroot + ynroot * znroot;
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

function complex_carlson_rf(x, y, z) {
    var xn = x;
    var yn = y;
    var zn = z;
    var a = cdivr(cadd(cadd(xn, yn), zn), 3);
    var eps = Math.max(cabs(csub(a, xn)), cabs(csub(a, yn)), cabs(csub(a, zn))) / errtol(x, y, z);
    var an = a;
    var power4 = 1;
    var n = 0;

    while (power4 * eps > cabs(an)) {
        var xnroot = csqrt(xn);
        var ynroot = csqrt(yn);
        var znroot = csqrt(zn);
        var lam = cmul(xnroot, ynroot);
        lam = cadd(lam, cmul(xnroot, znroot));
        lam = cadd(lam, cmul(ynroot, znroot));
        power4 *= 1/4;
        xn = cdivr(cadd(xn, lam), 4);
        yn = cdivr(cadd(yn, lam), 4);
        zn = cdivr(cadd(zn, lam), 4);
        an = cdivr(cadd(an, lam), 4);
        ++n;
    }
    var xndev = cdiv(rcmul(power4, csub(a, x)), an);
    var yndev = cdiv(rcmul(power4, csub(a, y)), an);
    var zndev = rcmul(-1, cadd(xndev, yndev));
    var ee2 = csub(cmul(xndev, yndev), rcmul(6, cmul(zndev, zndev)));
    var ee3 = cmul(cmul(xndev, yndev), zndev);
    //var s = 1 - ee2/10 + ee3/14 + ee2*ee2/24 - 3/44*ee2*ee3;
    var s = cadd(cdivr(cmul(ee2, ee2), 24), rcmul(-3/44, cmul(ee2, ee3)));
    s = cadd(s, cadd(cdivr(ee2, 10), cdivr(ee3, 14)));
    s = cadd(s, {re: 1, im: 0});
    return cdiv(s, csqrt(an));
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
        return Math.acosh(1/w);
    }
    return inverse_jacobi_cn(w, 1/m) / Math.sqrt(m);
        
}

function jacobi_cd(u, m) {
    return jacobi_cn(u, m) / jacobi_dn(u, m);
}

function inverse_jacobi_cd(w, m) {
    return inverse_jacobi_sn(Math.sqrt(1 - w*w) / Math.sqrt(1 - m*w*w), m);
}

function inverse_jacobi_cdk(w, k) {
    return inverse_jacobi_cd(w, k*k);
}

function complex_inverse_jacobi_sn(u, m) {
    var arg1 = cmul(u,u);
    arg1 = rcsub(1, arg1);

    var arg2 = rcsub(1, rcmul(m, cmul(u,u)));
    console.log("arg1");
    console.log(arg1)
    console.log("arg2");
    console.log(arg2)
        
    return cmul(u, complex_carlson_rf(arg1, arg2,{re: 1, im:0}));
}
