var filterType = "butterworth";
var context;
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

function selectFilter(type) {
  filterType = type;
}

function designFilter() {
  var passBand = document.getElementById("passband").value;
  var stopBand = document.getElementById("stopband").value;
  var passdB = document.getElementById("passdB").value;
  var stopdB = document.getElementById("stopdB").value;

  if (filterType == "butterworth") {
    butterworthCore(passBand, stopBand, passdB, stopdB, context.sampleRate);
  }
}

function butterworthPoles(passband, stopband, passdB, stopdB, sampleRate) {
  var f0 = 2 * Math.PI * passband / sampleRate;
  var f1 = 2 * Math.PI * stopband / sampleRate;
  var d0 = Math.pow(10, passdB / 10);
  var d1 = Math.pow(10, stopdB / 10);

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

function analogTeX(N, wc, terms) {
  var formula = "\\begin{align*} H_a(s) = & \\, " + Math.pow(wc, N) + "\\\\\n";

  if ((N & 1) == 1) {
    formula += "& \\frac{1}{(s + " + terms[0] + ")} \\\\\n";
    console.log("s + " + wc);
  }

  for (var k = (N & 1); k < terms.length; ++k) {
    formula += "& \\times \\frac{1}{s^2 + ";
    formula += terms[k][1] + "\\,s + ";
    formula += terms[k][2] + "}";
    formula += "\\\\\n";
  }

  formula += "\n\\end{align*}\n";

  return formula;
}

function digitalTeX(N, wc, terms) {
  var dFormula = "\\begin{align*} H_d(z) = &" + Math.pow(wc, N) + " \\\\\n";;
  if ((N & 1) == 1) {
    var b0 = (1 / (wc + 2));
    var a1 = ((wc - 2) / (wc + 2));
    dFormula += " & \\times \\frac{" + b0 + "(1+z^{-1})}";
    dFormula += "{1-" + (-a1) + "z^{-1}}\\\\\n";
  }
  for (var k = (N & 1); k < terms.length; ++k) {
    var b = terms[k][1];
    var c = terms[k][2];;
    var a0 = c + 2 * b + 4;

    console.log((1 / a0) + "*(1+1/z)^2/(1" + ((2 * c - 8) / a0) + "/z+" + ((c - 2 * b + 4) /
      a0) + "/z^2");

    dFormula += "& \\times";

    dFormula += "\\frac{" + (1 / a0) + "(1+z^{-1})^2}";
    dFormula += "{1-" + ((8 - 2 * c) / a0) + "\\,z^{-1}+" + ((c - 2 * b + 4) / a0) +
      "\\,z^{-2}}";
    dFormula += "\\\\\n";
  }
  dFormula += "\\end{align*}";

  return dFormula;
}

function webAudioFormula(N, wc, terms) {
  var waFormula;

  if (hasNewBiquadFilter) {
    waFormula = "<p>New biquad filters implemented; using BiquadFilter.</p>";
  } else {
    waFormula = "<p>New biquad filters not implemented; using IIRFilter instead.</p>";
  }

  waFormula += "<pre>\n";
  var totalGain = Math.pow(wc, N);;

  if ((N & 1) == 1) {
    var b0 = (1 / (wc + 2));
    var a1 = ((wc - 2) / (wc + 2));

    waFormula += "f0 = context.createIIRFilter([";
    waFormula += "[" + b0 + ", " + b0 + "], ";
    waFormula += "[" + 1 + ", " + a1 + "]);\n\n";
  }

  for (var k = (N & 1); k < terms.length; ++k) {
    var b = terms[k][1];
    var c = terms[k][2];;
    var a0 = c + 2 * b + 4;
    var f = "f" + k;

    if (hasNewBiquadFilter) {
      waFormula += f + " = context.createBiquadFilter();\n";
      waFormula += f + '.type = "lowpass";\n';

      var a1 = ((8 - 2 * c) / a0);
      var a2 = ((c - 2 * b + 4) / a0);
      var alpha = (1 - a2) / (1 + a2);
      var w0 = Math.acos((1 + alpha) * a1 / 2);
      var f0 = w0 * sampleRate / 2 / Math.PI;
      var Q = Math.sin(w0) / 2 / alpha;

      waFormula += f + ".frequency.value = " + f0 + "; // Hz\n"
      waFormula += f + ".Q.value = " + (20 * Math.log10(Q)) + "; // dB\n";

      //var gain = (1-Math.cos(w0))/2*a0;
      var gain = 2 / (1 - Math.cos(w0)) * (1 + alpha) / a0;
      totalGain *= gain;
      //waFormula += "g" + k + " = context.createGain();\n";
      //waFormula += "g" + k + ".gain.value = " + gain + ";\n";
      //waFormula += "\n";

    } else {
      waFormula += f + " = context.createIIRFilter(\n";
      waFormula += "        [";
      waFormula += (1 / a0) + ", " + (2 / a0) + ", " + (1 / a0) + "],\n";
      waFormula += "        ";
      waFormula += "[1, " + (-(8 - 2 * c) / a0) + ", " + ((c - 2 * b + 4) / a0) + "]);\n";
    }
  }

  for (var k = 1; k < terms.length; ++k) {
    waFormula += "f" + (k - 1) + ".connect(f" + k + ");\n";
    filters[k - 1].connect(filters[k]);
  }

  waFormula += "g = context.createGain();\n";
  waFormula += "g.gain.value = " + totalGain + ";\n";
  waFormula += "f" + (terms.length - 1) + ".connect(g);\n\n";
  waFormula += "g.connect(context.destination);\n";
  waFormula += "</pre>\n";
  return waFormula;
}

function createFilterGraph(N, wc, terms) {
  filters = new Array(terms.length);

  var totalGain = Math.pow(wc, N);;

  for (var k = 0; k < terms.length; ++k) {
    if (hasNewBiquadFilter) {
      var a1 = ((8 - 2 * c) / a0);
      var a2 = ((c - 2 * b + 4) / a0);
      var alpha = (1 - a2) / (1 + a2);
      var w0 = Math.acos((1 + alpha) * a1 / 2);
      var f0 = w0 * sampleRate / 2 / Math.PI;
      var Q = Math.sin(w0) / 2 / alpha;

      //var gain = (1-Math.cos(w0))/2*a0;
      var gain = 2 / (1 - Math.cos(w0)) * (1 + alpha) / a0;
      totalGain *= gain;

      filters[k] = context.createBiquadFilter();
      filters[k].type = "lowpass";
      filters[k].frequency.value = f0;
      filters[k].Q.value = 20 * Math.log10(Q);
    } else {
      if (terms[k].length == 3) {
        var b = terms[k][1];
        var c = terms[k][2];
        var a0 = c + 2 * b + 4;

        filters[k] = context.createIIRFilter(
          [1 / a0, 2 / a0, 1 / a0], [1, -(8 - 2 * c) / a0, (c - 2 * b + 4) / a0]);
      } else {
        var b0 = (1 / (wc + 2));
        var a1 = ((wc - 2) / (wc + 2));

        filters[k] = context.createIIRFilter([b0, b0], [1, a1]);
      }
    }
  }
  gain = context.createGain();
  gain.gain.value = totalGain;
  filters[terms.length - 1].connect(gain);
  gain.connect(context.destination);
}

function plotAnalogResponse(N, wc) {
  var freq = new Float32Array(1000);
  var mag = new Float32Array(freq.length);
  var gain = Math.pow(wc, N);

  for (var k = 0; k < freq.length; ++k) {
    freq[k] = k * context.sampleRate / 2 / freq.length;
    mag[k] = gain;
  }

  var dataAnalog = [];

  var cutoff = wc / Math.PI * context.sampleRate / 2;

  for (var k = 0; k < freq.length; ++k) {
    var r = 10 * Math.log10(1 / (1 + Math.pow(freq[k] / cutoff, 2 * N)));
    dataAnalog.push([freq[k], r]);
  }

  $.plot($("#graph-analog"), [{
    data: dataAnalog,
    label: "Magnitude response"
  }]);
}

function plotDigitalResponse(N, wc) {
  var freq = new Float32Array(1000);
  var mag = new Float32Array(freq.length);
  var phase = new Float32Array(freq.length);

  var gain = Math.pow(wc, N);
  var dataDigital = [];

  for (var k = 0; k < freq.length; ++k) {
    freq[k] = k * context.sampleRate / 2 / freq.length;
    mag[k] = gain;
    phase[k] = 0;
  }

  for (var k = 0; k < filters.length; ++k) {
    var m = new Float32Array(freq.length);
    var p = new Float32Array(freq, length);
    filters[k].getFrequencyResponse(freq, m, p);

    for (var n = 0; n < freq.length; ++n) {
      mag[n] *= m[n];
      phase[n] += p[n];
    }
  }

  // Plot the response

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

function butterworthCore(passband, stopband, passdB, stopdB, sampleRate) {
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

  if (hasNewBiquadFilter || hasIIRFilter) {
    createFilterGraph(N, wc, filterTerms);

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

  plotDigitalResponse(N, wc);
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

function init() {
  checkBiquadFilterQ().then(function (flag) {
    hasNewBiquadFilter = flag;
    context = new AudioContext();
    hasIIRFilter = context["createIIRFilter"] ? true : false;
  })
}
