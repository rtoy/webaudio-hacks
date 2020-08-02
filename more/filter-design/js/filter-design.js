var context;
var offline;
var sampleRate = 48000;
var hasNewBiquadFilter;
var hasIIRFilter;
var filters;
var freq;
var mag;
var phase;

var filterType = 'lowpass';
var plotType = 'dB';
var gain;
var osc;
var modGain;
var mod;

function setPlotType(type) {
  plotType = type;
}

function setFilterType(type) {
  filterType = type;
  if (type == 'lowpass') {
    var b1 = document.getElementById('band1');
    var b2 = document.getElementById('band2');
    var s1 = document.getElementById('band1-db');
    var s2 = document.getElementById('band2-db');
    b1.innerHTML = 'Passband (Hz)';
    s1.innerHTML = 'Passband attentuation, dB';
    b2.innerHTML = 'Stopband (Hz)';
    s2.innerHTML = 'Stopband attenuation, dB';

    // Swap the attenuation values.
    var tmp = document.getElementById('band1-db-value').value;
    document.getElementById('band1-db-value').value =
        document.getElementById('band2-db-value').value;
    document.getElementById('band2-db-value').value = tmp;

  } else if (type == 'highpass') {
    var b1 = document.getElementById('band1');
    var b2 = document.getElementById('band2');
    var s1 = document.getElementById('band1-db');
    var s2 = document.getElementById('band2-db');
    b1.innerHTML = 'Stopband (Hz)';
    s1.innerHTML = 'Stopband attentuation, dB';
    b2.innerHTML = 'Passband (Hz)';
    s2.innerHTML = 'Passband attenuation, dB';

    // Swap the attenuation values.
    var tmp = document.getElementById('band1-db-value').value;
    document.getElementById('band1-db-value').value =
        document.getElementById('band2-db-value').value;
    document.getElementById('band2-db-value').value = tmp;
  }
}

function designFilter(filterImplType) {
  if (filterType == 'lowpass') {
    designLowpassFilter(filterImplType);
  } else if (filterType == 'highpass') {
    designHighpassFilter(filterImplType);
  } else {
    alert('Not yet implemented: ' + filterType);
  }
}

function designLowpassFilter(filterImplType) {
  sampleRate = document.getElementById('samplerate').value;
  var passBand = Number(document.getElementById('band1-value').value);
  var stopBand = Number(document.getElementById('band2-value').value);
  var passdB = Number(document.getElementById('band1-db-value').value);
  var stopdB = Number(document.getElementById('band2-db-value').value);

  if (passBand <= 0 || stopBand <= 0 || passBand >= stopBand) {
    alert('Invalid passband or stopband frequencies');
    return;
  }

  var analogFilter =
      analogLowpassFilter(passBand, stopBand, passdB, stopdB, filterImplType);

  var aFormula = analogTeX(analogFilter);
  console.log(aFormula);

  var description = 'Analog lowpass ' + filterImplType + ' design of order ' +
      analogFilter.order;
  document.getElementById('analog-type').innerHTML = description;

  // var math = MathJax.Hub.getAllJax("analog-eq");
  // MathJax.Hub.Queue(["Text", math[0], aFormula]);
  MathJax.typesetPromise()
      .then(() => {
        const eqn = document.querySelector('#analog-eq');
        eqn.innerHTML = '$$' + aFormula + '$$';
        MathJax.typesetPromise();
      })
      .catch((err) => console.log(err.message));


  plotAnalogResponse(analogFilter);

  if (sampleRate <= 0) {
    alert('Sample rate must be positive!');
    return;
  }

  if (passBand >= sampleRate / 2 || stopBand >= sampleRate / 2) {
    alert(
        'Invalid pass band or stop band frequency exceeds Nyquist frequency.');
    return;
  }

  var digitalFilter = digitalLowpassFilter(
      passBand, stopBand, passdB, stopdB, sampleRate, filterImplType);
  var digitalTeXFormula = digitalTeX(digitalFilter);
  console.log(digitalTeXFormula);

  description = 'Digital lowpass ' + filterImplType + ' design of order ' +
      digitalFilter.order;
  description += ', sample rate ' + sampleRate + ' Hz';
  document.getElementById('digital-type').innerHTML = description;
  // math = MathJax.Hub.getAllJax("digital-eq");
  // MathJax.Hub.Queue(["Text", math[0], digitalTeXFormula]);
  MathJax.typesetPromise()
      .then(() => {
        const eqn = document.querySelector('#digital-eq');
        eqn.innerHTML = '$$' + digitalTeXFormula + '$$';
        MathJax.typesetPromise();
      })
      .catch((err) => console.log(err.message));

  plotDigitalResponse(digitalFilter, sampleRate);

  var webaudio = webAudioFilter(digitalFilter, sampleRate, filterImplType);
  console.log(webaudio);
  displayWebAudio(webaudio, {
    order: digitalFilter.order,
    sampleRate: sampleRate,
    filterType: 'lowpass',
    filterImplType: filterImplType,
    passBand: passBand,
    stopBand: stopBand,
    passAttenuation: passdB,
    stopAttenuation: stopdB
  });

  plotWebAudioResponse(webaudio, sampleRate, filterImplType);
}

function designHighpassFilter(filterImplType) {
  sampleRate = document.getElementById('samplerate').value;
  var passBand = Number(document.getElementById('band1-value').value);
  var stopBand = Number(document.getElementById('band2-value').value);
  var passdB = Number(document.getElementById('band2-db-value').value);
  var stopdB = Number(document.getElementById('band1-db-value').value);

  var analogFilter =
      analogHighpassFilter(passBand, stopBand, passdB, stopdB, filterImplType);

  console.log('Highpass:  lowpass equivalent:');
  console.log(analogFilter);

  var aFormula = analogTeX(analogFilter);
  console.log(aFormula);

  var description = 'Analog highpass ' + filterImplType + ' design of order ' +
      analogFilter.order;
  document.getElementById('analog-type').innerHTML = description;
  // var math = MathJax.Hub.getAllJax("analog-eq");
  // MathJax.Hub.Queue(["Text", math[0], aFormula]);
  MathJax.typesetPromise()
      .then(() => {
        const eqn = document.querySelector('#analog-eq');
        eqn.innerHTML = '$$' + aFormula + '$$';
        MathJax.typesetPromise();
      })
      .catch((err) => console.log(err.message));

  plotAnalogResponse(analogFilter);

  var digitalFilter = digitalHighpassFilter(
      passBand, stopBand, passdB, stopdB, sampleRate, filterImplType);
  var digitalTeXFormula = digitalTeX(digitalFilter);
  console.log(digitalTeXFormula);

  description = 'Digital highpass ' + filterImplType + ' design of order ' +
      digitalFilter.order;
  description += ', sample rate ' + sampleRate + ' Hz';
  document.getElementById('digital-type').innerHTML = description;
  // math = MathJax.Hub.getAllJax("digital-eq");
  // MathJax.Hub.Queue(["Text", math[0], digitalTeXFormula]);
  MathJax.typesetPromise()
      .then(() => {
        const eqn = document.querySelector('#analog-eq');
        eqn.innerHTML = '$$' + digitalTeXFormula + '$$';
        MathJax.typesetPromise();
      })
      .catch((err) => console.log(err.message));

  plotDigitalResponse(digitalFilter, sampleRate);

  var webaudio = webAudioFilter(digitalFilter, sampleRate, filterImplType);
  console.log(webaudio);
  displayWebAudio(webaudio, {
    order: digitalFilter.order,
    sampleRate: sampleRate,
    filterType: 'highpass',
    filterImplType: filterImplType,
    passBand: passBand,
    stopBand: stopBand,
    passAttenuation: passdB,
    stopAttenuation: stopdB
  });

  plotWebAudioResponse(webaudio, sampleRate, filterImplType);
}


function displayWebAudio(webaudioDesc, options) {
  // Generate code that implements the filter structure described by
  // webaudioDesc.
  var text = '<pre>\n';
  text += '<code class="javascript">\n';
  text += '// WebAudio ' + options.filterType + ' ' + options.filterImplType +
      ' design\n';
  text += '// Order = ' + options.order + '\n';
  text += '// Sample rate = ' + options.sampleRate + ' Hz\n';
  text += '// Passband = ' + options.passBand + ' Hz\n';
  text += '// Stopband = ' + options.stopBand + ' Hz\n';
  text += '// Passband attenuation = ' + options.passAttenuation + ' dB\n';
  text += '// Stopband attenuation = ' + options.stopAttenuation + ' dB\n';
  text += '\n';

  var filters = webaudioDesc.desc;
  for (var k = 0; k < filters.length; ++k) {
    var v = 'f' + k;
    text += v + ' = ';
    if (filters[k].filterType === 'biquad' && hasNewBiquadFilter) {
      text += 'context.createBiquadFilter();\n';
      text += v + '.type = "' + filters[k].biquadType + '";\n';
      text += v + '.frequency.value = ' + filters[k].f0 + ';\n';
      text += v + '.Q.value = ' + filters[k].Q + ';\n';
    } else {
      var g = filters[k].filterGain;
      text += 'context.createIIRFilter(\n';
      text += '       [';
      text += filters[k].top.map(x => x * g);
      text += '],\n';
      text += '       [';
      text += filters[k].bot;
      text += ']);\n';
    }
  }
  text += '\n';

  // Connect all the filters together
  for (var k = 1; k < filters.length; ++k) {
    text += 'f' + (k - 1);
    text += '.connect(' +
        'f' + k + ');\n';
  }

  // Create the gain term, if needed, and connect it
  if (webaudioDesc.totalGain != 1) {
    text += '\n';
    text += 'g = context.createGain();\n';
    text += 'g.gain.value = ' + webaudioDesc.totalGain + ';\n';
    text += 'f' + (filters.length - 1);
    text += '.connect(g);\n';
    text += 'g.connect(context.destination);\n';
  }

  text += '</code>\n';
  text += '</pre>\n';
  var element = document.getElementById('webaudio-eq');
  element.innerHTML = text;
  // hljs.highlightBlock(element);
}

function createGraph(webaudioDesc, Fs, filterType) {
  offline = new OfflineAudioContext(1, 1, Fs);
  var f = webaudioDesc.desc;
  filters = new Array(f.length);
  for (var k = 0; k < f.length; ++k) {
    if (f[k].filterType === 'biquad' && hasNewBiquadFilter) {
      filters[k] = offline.createBiquadFilter();
      filters[k].type = f[k].biquadType;
      filters[k].frequency.value = f[k].f0;
      filters[k].Q.value = f[k].Q;
    } else {
      var g = f[k].filterGain;
      filters[k] = offline.createIIRFilter(f[k].top.map(x => x * g), f[k].bot);
    }
  }

  for (var k = 1; k < f.length; ++k) {
    filters[k - 1].connect(filters[k]);
  }

  gain = offline.createGain();
  gain.gain.value = webaudioDesc.totalGain;
  filters[f.length - 1].connect(gain);
  gain.connect(offline.destination);
}

function plotWebAudioResponse(webaudioDesc, Fs, filterType) {
  try {
    context = new OfflineAudioContext(1, 1, Fs);
  } catch (e) {
    var warning =
        '<p>Could not create offline audio context with sample rate ' + Fs +
        ' Hz.</p>';
    warning += '<p>Thus the filter response cannot be plotted.</p>';
    document.getElementById('graph-webaudio').innerHTML = warning;
    return;
  };

  try {
    createGraph(webaudioDesc, Fs, filterType);
  } catch (e) {
    // If we get here, createGraph couldn't create the graph
    // because this implementation doesn't support IIRFilterNode.
    // Display a message where the graph would be.
    var warning =
        'This browser does not support IIRFilterNodes so no graph can be shown.';
    document.getElementById('graph-webaudio').innerHTML = warning;
    return;
  }

  var freq = new Float32Array(1000);

  for (var k = 0; k < freq.length; ++k) {
    freq[k] = k * sampleRate / 2 / freq.length;
  }

  var totalMag = new Float32Array(freq.length);
  var totalPhase = new Float32Array(freq.length);
  var mag = new Float32Array(freq.length);
  var phase = new Float32Array(freq.length);

  if (gain && gain.gain.value != 1) {
    totalMag.fill(gain.gain.value);
  } else {
    totalMag.fill(1);
  }
  for (var k = 0; k < filters.length; ++k) {
    filters[k].getFrequencyResponse(freq, mag, phase);
    for (var m = 0; m < mag.length; ++m) {
      totalMag[m] *= mag[m];
      totalPhase[m] += phase[m];
    }
  }

  var dataMag = [];
  var dataPhase = [];
  for (var k = 0; k < totalMag.length; ++k) {
    var r = plotType === 'dB' ? 20 * Math.log10(totalMag[k]) : totalMag[k];
    dataMag.push([freq[k], r]);
    dataPhase.push([freq[k], totalPhase[k] * 180 / Math.PI]);
  }

  let legendContainer = document.getElementById('graph-webaudio-legend');
  if (plotType == 'dB') {
    plotOptions = {
      grid: {hoverable: true},
      legend: {
        // position: 'se',
        show: true,
        container: legendContainer,
      },
      xaxis: {
        // mode: 'log',
        axisLabel: 'Freq (Hz)'
      },
      yaxes: [
        {min: -100, max: 5, axisLabel: 'Mag'},
        {alignTicksWithAxis: 1, position: 'right', axisLabel: 'Phase (deg)'}
      ]
    };

  } else {
    plotOptions = {
      grid: {hoverable: true},
      legend: {
        // position: 'se',
        show: true,
        container: legendContainer,
      },
      xaxis: {
          // mode: 'log',
      },
      yaxes: [
        {}, {
          position: 'right',
          alignTicksWithAxis: 1,
        }
      ]
    };
  }

  $.plot(
      $('#graph-webaudio'),
      [
        {
          data: dataMag,
          label: 'Magnitude response',
          lines: {lineWidth: 2},
          color: 'red'
        },
        {
          data: dataPhase,
          label: 'Phase response',
          lines: {lineWidt: 2},
          color: 'green',
          yaxis: 2
        }
      ],
      plotOptions);
}

function plotAnalogResponse(filter) {
  var freq = new Float32Array(1000);

  for (var k = 0; k < freq.length; ++k) {
    freq[k] = k * sampleRate / 2 / freq.length;
  }

  // var {mag, phase} = analogResponse(filter, freq);
  var response = analogResponse(filter, freq);
  mag = response.mag;
  phase = response.phase;

  console.log(freq);
  console.log(mag);

  var analogMag = [];
  var analogPhase = [];
  for (var k = 0; k < freq.length; ++k) {
    var r = plotType === 'dB' ? 20 * Math.log10(mag[k]) : mag[k];
    analogMag.push([freq[k], r]);
    analogPhase.push([freq[k], (phase[k] * 180 / Math.PI) % 180]);
  }

  var plotOptions;

  let legendContainer = document.getElementById('graph-analog-legend');
  if (plotType == 'dB') {
    plotOptions = {
      grid: {hoverable: true},
      legend: {
        // position: 'se',
        show: true,
        container: legendContainer,
      },
      xaxis: {
        // mode: 'log',
        axisLabel: 'Freq (Hz)'
      },
      yaxes: [
        {min: -80, axisLabel: 'Mag'},
        {position: 'right', alignTicksWithAxis: 1, axisLabel: 'Phase (deg)'}
      ]
    };
  } else {
    plotOptions = {yaxes: [{}, {position: 'right'}]};
  }

  $.plot(
      $('#graph-analog'),
      [
        {
          data: analogMag,
          label: 'Magnitude response',
          lines: {lineWidth: 2},
          color: 'red'
        },
        {
          data: analogPhase,
          label: 'Phase response',
          label: 'Magnitude response',
          lines: {lineWidth: 2},
          color: 'green',
          yaxis: 2
        }
      ],
      plotOptions);
}

function plotDigitalResponse(filter, Fs) {
  var freq = new Float32Array(1000);

  for (var k = 0; k < freq.length; ++k) {
    freq[k] = k * sampleRate / 2 / freq.length;
  }

  // var {mag, phase} = digitalResponse(filter, freq, Fs);
  var response = digitalResponse(filter, freq, Fs);
  mag = response.mag;
  phase = response.phase;

  // Plot the response
  var digitalMag = [];
  var digitalPhase = [];
  for (var k = 0; k < freq.length; ++k) {
    var r = plotType === 'dB' ? 20 * Math.log10(mag[k]) : mag[k];
    digitalMag.push([freq[k], r]);
    digitalPhase.push([freq[k], phase[k] * 180 / Math.PI]);
  }

  var plotOptions;

  let legendContainer = document.getElementById('graph-digital-legend');
  if (plotType == 'dB') {
    plotOptions = {
      legend: {
        // position: 'se',
        show: true,
        container: legendContainer,
      },
      xaxis: {axisLabel: 'Freq (Hz)'},

      yaxes: [
        {min: -80, axisLabel: 'Mag'},
        {alignTicksWithAxis: 1, position: 'right', axisLabel: 'Phase (deg)'}
      ]
    };

  } else {
    plotOptions = {
      grid: {hoverable: true},
      legend: {
        // position: 'ne',
        show: true,
        // container: legendContainer,
      },
      yaxes: [
        {}, {
          position: 'right',
          alignTicksWithAxis: 1,
        }
      ]
    };
  }

  $.plot(
      $('#graph-digital'),
      [
        {
          data: digitalMag,
          label: 'Magnitude response',
          lines: {lineWidth: 2},
          color: 'red'
        },
        {
          data: digitalPhase,
          label: 'Phase response',
          lines: {lineWidt: 2},
          color: 'green',
          yaxis: 2
        }
      ],
      plotOptions);
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

  return context.startRendering().then(function(buffer) {
    return Math.max(...buffer.getChannelData(0)) !== 0;
  });
}

// Initialization.  Determine if this browser supports the biquad
// filter with the new Q interpretation and if the browser supports
// IIRFilters.
function init() {
  checkBiquadFilterQ().then(function(flag) {
    hasNewBiquadFilter = flag;
    context = new AudioContext();
    hasIIRFilter = context['createIIRFilter'] ? true : false;
  })
}
