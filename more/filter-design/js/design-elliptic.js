function designElliptic(fp, fs, Gp, Gs, sampleRate) {
    var wp = 2 * Math.PI * fp;
    var ws = 2 * Math.PI * fs;
    var ep = Math.sqrt(1/(Gp*Gp) - 1);
    var es = Math.sqrt(1/(Gs*Gs) - 1);

    console.log("ep = " + ep + ", es = " + es);

    var k = wp/ws;
    var k1 = ep/es;

    console.log("k = " + k + ", k1 = " + k1);

    var K = elliptic_kc(k*k);
    var Kp = elliptic_kc(1-k*k);
    console.log("K = " + K + ", K' = " + Kp);

    var K1 = elliptic_kc(k1*k1);
    var K1p = elliptic_kc(1-k1*k1);
    console.log("K1 = " + K1 + ", K1' = " + K1p);

    var N = (K1p/K1)/(Kp/K);
    console.log("N = " + N + ", " + Math.ceil(N));
    N = Math.ceil(N);

    k = ellipticDeg(N, K1, K1p);
    K = elliptic_kc(k*k);
    console.log("new k = " + k + ", K = " + K);

    var fs_new = fp / k;
    console.log("fs_new = " + fs_new);
    var L = Math.floor(N/2);
    var r = N - 2*L;

    var zeta_i = new Array(L);
    for (var n = 1; n <= L; ++n) {
        zeta_i[n-1] = jacobi_cd((2*n-1)/N*K, k*k);
    }
    console.log("zeta_i = ");
    console.log(zeta_i);

    var za = zeta_i.map(function (z) {return {im: wp / (k*z)}; });
    console.log("za = ");
    console.log(za);

    var v0 = inverse_jacobi_sni(1 / ep, k1*k1);
    v0 = -v0.im / N / K1;
    console.log("v0 = ");
    console.log(v0);
    
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
