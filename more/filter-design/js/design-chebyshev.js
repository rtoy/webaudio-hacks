function designChebyshevFilter(passBand, stopBand, passdB, stopdB, sampleRate) {
    var result = chebyshevPoles(passBand, stopBand, passdB, stopdB, context.sampleRate);
    var f0 = result.f0;
    var N = result.N;
    var eps = result.eps;
    var angle = result.poleAngles;

    var filterTerms = computeChebyshevFilter(N, f0, eps, angle);
    console.log("filterTerms = ");
    console.log(filterTerms);

    var analogTeXFormula = analogTeX(N, Math.pow(filterTerms[0], 1 / N), filterTerms[1]);
    var math = MathJax.Hub.getAllJax("analog-eq");
    MathJax.Hub.Queue(["Text", math[0], N]);
    MathJax.Hub.Queue(["Text", math[1], analogTeXFormula]);
    MathJax.Hub.Queue(["Text", math[2], passBand]);

    plotAnalogChebyshev(N, passBand | 0, eps, filterTerms);

    var digitalTeXFormula = digitalTeX(N, Math.pow(filterTerms[0], 1 / N), filterTerms[1]);
    math = MathJax.Hub.getAllJax("digital-eq");
    MathJax.Hub.Queue(["Text", math[0], context.sampleRate]);
    MathJax.Hub.Queue(["Text", math[1], digitalTeXFormula]);

    if (hasNewBiquadFilter || hasIIRFilter) {
        createFilterGraph(N, Math.pow(filterTerms[0], 1 / N), filterTerms[1]);
    }

    plotDigitalResponse(N, Math.pow(filterTerms[0], 1 / N));

    var webaudioFormula = webAudioFormula(N, Math.pow(filterTerms[0], 1 / N), filterTerms[1]);
    document.getElementById("webaudio").innerHTML = webaudioFormula;
}

function chebyshev(n, x) {
    if (Math.abs(x) <= 1)
        return Math.cos(n * Math.acos(x));
    return Math.cosh(n * Math.acosh(x));
}

function chebyshevPoles(passBand, stopBand, passdB, stopdB, sampleRate) {
    var f0 = 2 * Math.PI * passBand / sampleRate;
    var f1 = 2 * Math.PI * stopBand / sampleRate;
    var d0 = Math.pow(10, passdB / 10);
    var d1 = Math.pow(10, stopdB / 10);

    // 1/(1+e^2) = 1/d0;
    // e = sqrt(d0-1)
    var eps = Math.sqrt(d0 - 1);

    // 1/(1+e^2*V(n,x)^2) <= 1/d1.
    var N = Math.acosh(Math.sqrt(d1 - 1) / eps) / Math.acosh(2 * Math.tan(f1 / 2) / (2 * Math.tan(
        f0 / 2)));

    console.log("e = " + eps + ", N = " + N);

    N = Math.ceil(N);

    // Find the angles for the roots.
    // For N even, the angles are (k+1/2)/N*pi, k = 0, 1, ... N/2-1
    // For N odd, the angles are k/N*%pi, k = 0, 1, ..., (N-1)/2

    var angle = Array(Math.floor(N / 2));

    if ((N & 1) == 0) {
        for (var k = 0; k <= N / 2 - 1; ++k) {
            angle[k] = (k + 1 / 2) / N * Math.PI;
        }
    } else {
        for (var k = 0; k <= (N - 1) / 2; ++k) {
            angle[k] = k / N * Math.PI;
        }
    }

    return {
        N: N,
        f0: f0,
        eps: eps,
        poleAngles: angle
    };
}

function computeChebyshevFilter(N, f0, eps, angle) {
    var terms = new Array(Math.floor(N / 2));
    var a = Math.sinh(Math.asinh(1 / eps) / N);
    var b = Math.cosh(Math.asinh(1 / eps) / N);

    var numer = ((N & 1) == 1) ? 1 : 1 / Math.sqrt(1 + eps * eps);;
    console.log("angles");
    console.log(angle);
    for (var k = 0; k < angle.length; ++k) {
        // Poles are located at
        // x = -a*f0*cos(angle)
        // y = b*f0*sin(angle)
        var x = -a * f0 * Math.cos(angle[k]);
        if ((N & 1) == 1 && (k == 0)) {
            // For an odd filter, the first term is linear.
            terms[k] = [1, -x];
            numer *= -x;
        } else {
            var y = b * f0 * Math.sin(angle[k]);
            var r2 = x * x + y * y;
            console.log("pole " + k + ": (" + x + ", " + y + ")");
            terms[k] = [1, -2 * x, r2];

            numer *= r2;
        }
    }

    console.log("numerator: " + numer);
    return [numer, terms];
}

function plotAnalogChebyshev(N, passBand, eps, filterTerms) {
    var freq = new Float32Array(1000);
    var mag = new Float32Array(freq.length);
    var gain = filterTerms[0];

    for (var k = 0; k < freq.length; ++k) {
        freq[k] = k * context.sampleRate / 2 / freq.length;
        mag[k] = gain;
    }

    var dataAnalog = [];

    for (var k = 0; k < freq.length; ++k) {
        var r = 1 / (1 + eps * eps * Math.pow(chebyshev(N, freq[k] / passBand), 2));
        dataAnalog.push([freq[k], 10 * Math.log10(r)]);
    }

    $.plot($("#graph-analog"), [{
        data: dataAnalog,
        label: "Magnitude response"
    }]);
}
