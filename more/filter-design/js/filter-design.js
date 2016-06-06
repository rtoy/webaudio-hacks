var context;
var offline;
var sampleRate = 48000;
var hasNewBiquadFilter;
var hasIIRFilter;
var filters;
var freq;
var mag;
var phase;

var gain;
var osc;
var modGain;
var mod;

function designFilter(filterType) {
    sampleRate = document.getElementById("samplerate").value;
    var passBand = document.getElementById("passband").value;
    var stopBand = document.getElementById("stopband").value;
    var passdB = document.getElementById("passdB").value;
    var stopdB = document.getElementById("stopdB").value;

    var analogFilter = analogLowpassFilter(passBand, stopBand, passdB, stopdB, filterType);
    var digitalFilter = digitalLowpassFilter(passBand, stopBand, passdB, stopdB, sampleRate, filterType);

    var aFormula = analogTeX(analogFilter);
    console.log(aFormula);

    var description = "Analog " + filterType + " design of order " + analogFilter.order;
    document.getElementById("analog-type").innerHTML = description;
    var math = MathJax.Hub.getAllJax("analog-eq");
    MathJax.Hub.Queue(["Text", math[0], aFormula]);

    plotAnalogResponse(analogFilter);

    var digitalTeXFormula = digitalTeX(digitalFilter);
    console.log(digitalTeXFormula);

    description = "Digital " + filterType + " design of order " + digitalFilter.order;
    description += ", sample rate " + sampleRate + " Hz";
    document.getElementById("digital-type").innerHTML = description;
    math = MathJax.Hub.getAllJax("digital-eq");
    MathJax.Hub.Queue(["Text", math[0], digitalTeXFormula]);

    plotDigitalResponse(digitalFilter, sampleRate);

    var webaudio = webAudioFilter(digitalFilter, sampleRate, filterType);
    console.log(webaudio);
    displayWebAudio(webaudio);

    plotWebAudioResponse(webaudio, sampleRate, filterType);
}

function displayWebAudio(webaudioDesc) {
    // Generate code that implements the filter structure described by webaudioDesc.
    var text = "<pre>\n";
    var filters = webaudioDesc.desc;
    for (var k = 0; k < filters.length; ++k) {
	var v = "f" + k;
	text += v + " = ";
	if (filters[k].filterType === "biquad" && hasNewBiquadFilter) {
	    text += "context.createBiquadFilter();\n";
	    text += v + ".type = " + filters[k].biquadType + ";\n";
	    text += v + ".frequency.value = " + filters[k].f0 + ";\n";
	    text += v + ".Q.value = " + filters[k].Q + ";\n";
	} else {
	    text += "context.createIIRFilter(\n";
	    text += "       [";
	    text += filters[k].top;
	    text += "],\n";
	    text += "       [";
	    text += filters[k].bot;
	    text += "]);\n";
	}
    }
    text += "\n";

    // Connect all the filters together
    for (var k = 1; k < filters.length; ++k) {
	text += "f" + (k - 1);
	text += ".connect(" + "f" + k + ");\n";
    }

    // Create the gain term and conenct it
    text += "\n";
    text += "g = context.createGain();\n";
    text += "g.gain.value = " + webaudioDesc.totalGain + ";\n";
    text += "f" + (filters.length - 1);
    text += ".connect(g);\n";
    text += "g.connect(context.destination);\n";
    text += "</pre>\n";

    document.getElementById("webaudio").innerHTML = text;
}

function createGraph(webaudioDesc, Fs, filterType) {
    offline = new OfflineAudioContext(1, 1, Fs);
    var f = webaudioDesc.desc;
    filters = new Array(f.length);
    for (var k = 0; k < f.length; ++k) {
	if (f[k].filterType === "biquad" && hasNewBiquadFilter) {
	    filters[k] = offline.createBiquadFilter();
	    filters[k].type = f[k].biquadType;
	    filters[k].frequency.value = f[k].f0;
	    filters[k].Q.value = f[k].Q;
	} else {
	    filters[k] = offline.createIIRFilter(f[k].top, f[k].bot);
	}
    }

    for (var k = 1; k < f.length; ++k) {
	filters[k-1].connect(filters[k]);
    }

    gain = offline.createGain();
    gain.gain.value = webaudioDesc.totalGain;
    filters[f.length - 1].connect(gain);
    gain.connect(offline.destination);
}

function plotWebAudioResponse(webaudioDesc, Fs, filterType) {
    createGraph(webaudioDesc, Fs, filterType);

    var freq = new Float32Array(1000);

    for (var k = 0; k < freq.length; ++k) {
        freq[k] = k * sampleRate / 2 / freq.length;
    }

    var totalMag = new Float32Array(freq.length);
    var mag = new Float32Array(freq.length);
    var phase = new Float32Array(freq.length);

    totalMag.fill(gain.gain.value);
    for (var k = 0; k < filters.length; ++k) {
	filters[k].getFrequencyResponse(freq, mag, phase);
	for (var m = 0; m < mag.length; ++m) {
	    totalMag[k] *= mag[k];
	}
    }

    var data = [];
    for (var k = 0; k < totalMag.length; ++k) {
	data.push([freq[k], 20*Math.log10(totalMag[k])]);
    }
    
    $.plot($("#graph-webaudio"), [{
        data: data,
        label: "Magnitude response"
    }], {
        yaxes: [{
            min: -80
        }]
    });

}

// Find the poles of a Butterworth filter, returning the order, the
// (normalized cutoff frequency), the angles of the poles.  All the
// poles of a Butterworth filter lie on a circle, so only the angle
// and the radius is needed.  The radius is given by the cutoff
// frequency.
function butterworthPoles(passband, stopband, passdB, stopdB, sampleRate) {
    var f0 = 2 * Math.PI * passband / sampleRate;
    var f1 = 2 * Math.PI * stopband / sampleRate;
    var d0 = Math.pow(10, passdB / 10);
    var d1 = Math.pow(10, stopdB / 10);

    // Since we're using a bilinear transformation to determine design
    // the digital filter, we need to prewarp the frequencies to
    // determine the order and cutoff frequency.
    var N = 0.5 * Math.log((d1 - 1) / (d0 - 1)) / Math.log(Math.tan(f1 / 2) / Math.tan(f0 / 2));
    console.log("N = ", N);

    N = Math.ceil(N);
    var wc = 2 * Math.tan(f1 / 2) / Math.pow(d1 - 1, 1 / (2 * N));
    console.log("N = ", N, "wc = ", wc);

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
        wc: wc,
        f0: f0,
        f1: f1,
        d0: d0,
        d1: d1,
        poleAngles: angle
    };
}

function computeButterworthFilter(N, wc, angle) {
    // From the angles for the poles, the actual poles of the filter are then given by
    //
    //    wc*(-cos(a) +/- i*sin(a))
    //
    // The quadratic for this is then
    //
    //   s^2 + 2*wc*cos(a)*s + wc^2
    //
    // However, when N is odd, the first pole is real at -wc, so instead of a quadratic we have
    // a linear term (s + wc).
    //
    // Return an array consisting of the formula for each factor.
    //
    var terms = new Array(angle.length);
    for (var k = 0; k < angle.length; ++k) {
        if (k == 0 && (N & 1) == 1) {
            terms[k] = [1, wc];
        } else {
            terms[k] = [1, 2 * wc * Math.cos(angle[k]), wc * wc]
        }
    }

    return terms;
}

function applyBilinearTransform(N, wc, terms) {
    var xfrm = new Array(terms.length);
    console.log("applyBilinear to ");
    console.log(terms);

    for (var k = 0; k < terms.length; ++k) {
        var term = terms[k];

        if (term.length == 2) {
            // Linear 
            var b0 = (1 / (wc + 2));
            var a1 = ((wc - 2) / (wc + 2));
            xfrm[k] = [b0, [1, a1]];
        } else {
            // Quadratic
            var b = term[1];
            var c = term[2];
            var a0 = c + 2 * b + 4;

            xfrm[k] = [1 / a0, [1, (2 * c - 8) / a0, (c - 2 * b + 4) / a0]];
        }
    }

    return xfrm;
}

function plotAnalogResponse(filter) {
    var freq = new Float32Array(1000);

    for (var k = 0; k < freq.length; ++k) {
        freq[k] = k * sampleRate / 2 / freq.length;
    }

    var mag = analogResponse(filter, freq);

    console.log(freq);
    console.log(mag);

    var dataAnalog = [];
    for (var k = 0; k < freq.length; ++k) {
        var r = 20 * Math.log10(mag[k]);
        dataAnalog.push([freq[k], r]);
    }

    $.plot($("#graph-analog"), [{
        data: dataAnalog,
        label: "Magnitude response"
    }], {
        yaxes: [{
            min: -80
        }]
    });
}

function plotDigitalResponse(filter, Fs) {
    var freq = new Float32Array(1000);

    for (var k = 0; k < freq.length; ++k) {
        freq[k] = k * sampleRate / 2 / freq.length;
    }

    var mag = digitalResponse(filter, freq, Fs);

    // Plot the response
    var dataDigital = [];
    for (var k = 0; k < freq.length; ++k) {
        dataDigital.push([freq[k], 20 * Math.log10(mag[k])]);
    }

    $.plot($("#graph-digital"), [{
        data: dataDigital,
        label: "Magnitude response"
    }], {
        yaxes: [{
            min: -80
        }]
    });
}

function designButterworthFilter(passband, stopband, passdB, stopdB, sampleRate) {
    var result = butterworthPoles(passband, stopband, passdB, stopdB, sampleRate);
    var f0 = result.f0;
    var f1 = result.f1;
    var d0 = result.d0;
    var d1 = result.d1;
    var angle = result.poleAngles;
    var wc = result.wc;
    var N = result.N;

    var filterTerms = computeButterworthFilter(N, wc, angle);
    var analogTeXFormula = analogTeX(N, wc, filterTerms);

    console.log(analogTeXFormula);

    var math = MathJax.Hub.getAllJax("analog-eq");
    MathJax.Hub.Queue(["Text", math[0], N]);
    MathJax.Hub.Queue(["Text", math[1], analogTeXFormula]);
    MathJax.Hub.Queue(["Text", math[2], wc / Math.PI * sampleRate / 2]);

    // Convert the analog filter to digital filter using a bilinear transform.

    var digitalTeXFormula = digitalTeX(N, wc, filterTerms);
    console.log(digitalTeXFormula);

    math = MathJax.Hub.getAllJax("digital-eq");
    MathJax.Hub.Queue(["Text", math[0], context.sampleRate]);
    MathJax.Hub.Queue(["Text", math[1], digitalTeXFormula]);

    var totalGain = 1;
    if (hasNewBiquadFilter || hasIIRFilter) {
        totalGain = createFilterGraph(N, Math.pow(wc, N), filterTerms);

        osc = context.createOscillator();
        osc.type = "sine";
        osc.frequency.value = context.sampleRate / 4;;
        osc.connect(filters[0]);

        mod = context.createOscillator();
        mod.type = "sine";
        mod.frequency.value = 100;
        modGain = context.createGain();
        modGain.gain.value = context.sampleRate / 4;
        mod.connect(modGain);
        modGain.connect(osc.detune);
        if (0) {
            mod.start();
            osc.start();
        }
    }
    var webaudioFormula = webAudioFormula(N, wc, filterTerms);
    document.getElementById("webaudio").innerHTML = webaudioFormula;

    plotAnalogResponse(N, wc);

    plotDigitalResponse(N, totalGain);
}

// Check the Q implementation and return a promise.
function checkBiquadFilterQ() {
    'use strict';
    var context = new OfflineAudioContext(1, 128, 48000);
    var osc = context.createOscillator();
    var filter1 = context.createBiquadFilter();
    var filter2 = context.createBiquadFilter();
    var inverter = context.createGain();
    osc.type = 'sawtooth';
    osc.frequency.value = 8 * 440;
    inverter.gain.value = -1;
    // each filter should get a different Q value
    filter1.Q.value = -1;
    filter2.Q.value = -20;
    osc.connect(filter1);
    osc.connect(filter2);
    filter1.connect(context.destination);
    filter2.connect(inverter);
    inverter.connect(context.destination);
    osc.start();

    return context.startRendering().then(function (buffer) {
        return Math.max(...buffer.getChannelData(0)) !== 0;
    });
}

// Initialization.  Determine if this browser supports the biquad
// filter with the new Q interpretation and if the browser supports
// IIRFilters.
function init() {
    checkBiquadFilterQ().then(function (flag) {
        hasNewBiquadFilter = flag;
        context = new AudioContext();
        hasIIRFilter = context["createIIRFilter"] ? true : false;
    })
}
