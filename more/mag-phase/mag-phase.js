let context;
let filter;
let outputGain;
let sound;
let cutoff = 350;
let q = 10.0;    // in dB
let gain = 0.0;  // in dB
let source = null;

// Lowest frequency we want to plot.  This should be a power of 10 for
// the nicest graph.
const lowestFrequency = 10;

let sampleRate = 44100.0;
let nyquist = 0.5 * sampleRate;

let plot;

function dBFormatter(v, axis) {
  return v.toFixed(axis.tickDecimals) + ' dB';
}

function degFormatter(v, axis) {
  return v.toFixed(axis.tickDecimals) + ' deg';
}

function toms463(xmin, xmax, n) {
  // From TOMS 463
  var sqr = [1.414214, 3.162278, 7.071068];
  var vint = [1, 2, 5, 10];
  var del = 0.00002
  // Arbitrarily use 4 intervals.
  var a = (xmax - xmin) / 4;
  var al = Math.log(a) / Math.LN10;
  var nal = Math.floor(al);
  if (a < 1) {
    nal = nal - 1;
  }
  var b = a / Math.pow(10, nal);
  var i = 4;
  for (i = 0; i < 3; ++i) {
    if (b < sqr[i]) {
      break;
    }
  }
  var dist = vint[i] * Math.pow(10, nal);
  var fm1 = xmin / dist;
  var m1 = Math.floor(fm1);
  if (fm1 < 0) {
    m1 = m1 - 1;
  }
  if (Math.abs(m1 + 1 - fm1) < del) {
    m1 = m1 + 1;
  }
  var xminp = dist * m1;
  var fm2 = xmax / dist;
  var m2 = Math.floor(fm2 + 1);
  if (fm2 < -1) {
    m2 = m2 - 1;
  }
  if (Math.abs(fm2 + 1 - m2) < del) {
    m2 = m2 - 1;
  }
  var xmaxp = dist * m2;
  if (xminp > xmin) {
    xminp = xmin;
  }
  if (xmaxp < xmax) {
    xmaxp = xmax;
  }
  return [xminp, xmaxp, dist];
}
function tickScale(axis) {
  // Compute scale
  var tickInfo = toms463(axis.min, axis.max, 4);

  // Generate ticks now.
  var ticks = [];
  var val = tickInfo[0];
  while (val <= tickInfo[1]) {
    ticks.push(val);
    val = val + tickInfo[2];
  }
  return ticks;
}

function drawCurve() {
  // Lowest frequency
  var width = 1000;

  var freq = new Float32Array(width);
  var magResponse = new Float32Array(width);
  var phaseResponse = new Float32Array(width);

  // Logarithmically sample between lowest frequency and Nyquist by
  // uniformly sampling between the logs of the frequencies.

  let delta = Math.log10(nyquist / lowestFrequency) / width;
  let logLowest = Math.log10(lowestFrequency);

  for (var k = 0; k < width; ++k) {
    var f = logLowest + k * delta;
    freq[k] = Math.pow(10, f);
  }

  filter.getFrequencyResponse(freq, magResponse, phaseResponse);

  var magData = [];
  var phaseData = [];

  for (var k = 0; k < width; ++k) {
    db = 20.0 * Math.log(magResponse[k]) / Math.LN10;
    phaseDeg = 180 / Math.PI * phaseResponse[k];
    magData.push([freq[k], db]);
    phaseData.push([freq[k], phaseDeg]);
  }

  // Figure out the y axis range based on the filter type.

  var type = filter.type;

  switch (type) {
    case 'lowpass':
      magmin = -80;
      magmax = 40;
      phasemin = -200;
      phasemax = 0;
      break;
    case 'highpass':
      magmin = -80;
      magmax = 20;
      phasemin = 0;
      phasemax = 180;
      break;
    case 'bandpass':
      magmin = -80;
      magmax = 0;
      phasemin = -180;
      phasemax = 180;
      break;
    case 'lowshelf': {
      // Get the limits from the gain slider
      let slider = document.getElementById('gainSlider');
      magmin = slider.min;
      magmax = slider.max;
      phasemin = -180;
      phasemax = 180;
    } break;
    case 'highshelf': {
      // Get the limits from the gain slider
      let slider = document.getElementById('gainSlider');
      magmin = slider.min;
      magmax = slider.max;
      phasemin = -180;
      phasemax = 180;
    } break;
    case 'peaking': {
      // Get the limits from the gain slider
      let slider = document.getElementById('gainSlider');
      magmin = slider.min;
      magmax = slider.max;
      phasemin = -90;
      phasemax = 90;
    } break;
    case 'notch':
      magmin = -60;
      magmax = 0;
      phasemin = -100;
      phasemax = 100;
      break;
    case 'allpass':
      magmin = -1;
      magmax = 1;
      phasemin = -180;
      phasemax = 180;
      break;
    default:
      console.log(`Unknown filter type: ${type}`)
      break;
  }

  console.log('magmin: ' + magmin + ' magmax: ' + magmax);
  plot = $.plot(
      $('#graph'),
      [
        {data: magData, label: 'Mag (dB)', lines: {lineWidth: 3}, color: 'red'},
        {
          data: phaseData,
          label: 'Phase (deg)',
          lines: {lineWidth: 3},
          color: 'green',
          yaxis: 2
        }
      ],
      {
        xaxis: {
          mode: 'log',
          ticks: 10,
          showMinorTicks: true,
          axisLabel: 'Freq (Hz)'
        },
        // xaxes: [ { ticks : tickScale } ],
        yaxes: [
          {
            position: 'left',
            axisLabel: 'Mag (dB)',
            // tickFormatter: dBFormatter,
            min: magmin,
            max: magmax,
            // ticks : tickScale
            autoScale: 'none',
            tickDecimals: 1,
          },
          {
            position: 'right',
            axisLabel: 'Phase (deg)',
            // align if we are to the right
            alignTicksWithAxis: 1,
            // tickFormatter: degFormatter,
            min: phasemin,
            max: phasemax,
            // ticks : tickScale
            autoScale: 'none'
          }
        ],
        legend: {position: 'ne', show: true},
        grid: {hoverable: true}
      });
  $('#graph').bind('plothover', (event, pos, item) => {
    if (!pos.x || !pos.y) {
      return;
    }
    let str = `(${pos.x.toFixed(2)}, ${pos.y.toFixed(2)})`;
    $('#hoverdata').text(str);
    if (item) {
      let x = item.datapoint[0].toFixed(2);
      let y = item.datapoint[1].toFixed(2);
      $('#tooltip')
          .html(`${item.series.label} at ${x} = ${y}`)
          .css({top: item.pageY + 5, left: item.pageX + 5})
          .fadeIn(200);
    } else {
      $('#tooltip').stop().hide();
    }
  });
  $('#graph').bind('plothovercleanup', (event, pos, item) => {
    $('#tooltip').hide();
  });

  x_axis = plot.getAxes().xaxis;
  // Add a div to the x axis so we can figure some things about the axis
  let box = x_axis.box;

  $('<div class=\'axisTarget\' style=\'position:absolute; left:' + box.left +
    'px; top:' + box.top + 'px; width:' + box.width +
    'px; height:' + box.height + 'px\'></div>')
      .data('axis.direction', x_axis.direction)
      .data('axis.n', x_axis.n)
      .css({backgroundColor: '#f00', opacity: 0, cursor: 'pointer'})
      .appendTo(plot.getPlaceholder());

  adjustSliderPositions();
}

function adjustSliderPositions() {
  console.log('adjustSliderPosition');

  let box = x_axis.box;
  console.log(x_axis);

  console.log(box);

  // Now get the dimensions of the axisTarget and use that to adjust the size of
  // the sliders.
  let axisTarget = document.querySelector('.axisTarget')
  axisTarget.style.left = `${box.left}px`;
  axisTarget.style.top = `${box.top}px`;
  axisTarget.style.width = `${box.width}px`;
  axisTarget.style.height = `${box.height}px`;

  // Update controls container to move sliders to the right a bit.
  // But the axisTargetStyle.left value is just a little too far to
  // left. So add a fudge factor.  Haven't figured out how to get the
  // actual position of the left axis.
  let controls = document.getElementById('controls');
  const leftExtra = 10;

  controls.style.paddingLeft = `${box.left + leftExtra}px`;


  // Set the slider widths appropriately
  let sliders = document.getElementsByClassName('slider-bar');

  Array.prototype.forEach.call(sliders, (slider) => {
    slider.style.width = `${box.width - 30}px`;
  });
}

function stopSound() {
  if (source) {
    source.stop(0);
    source = null;
  }
}

function loadSound(url) {
  // Load asynchronously
  var request = new XMLHttpRequest();
  request.open('GET', url, true);
  request.responseType = 'arraybuffer';
  request.onload = function() {
    context.decodeAudioData(
        request.response,
        function(buffer) {
          sound = buffer;
          if (source) {
            source.stop(0);
            source = null;
          }
          context.resume().then(() => {
            source = context.createBufferSource();
            source.connect(filter);
            source.loop = true;
            source.buffer = buffer;
            source.start(0);
          });
        },
        function() {
          console.log('error decoding file.')
        });
  };


  request.send();
}

function cutoffHandler(event, ui) {
  console.log('cutoffHandler ' + event + ' ' + ui.value);
  let cutoff = Math.pow(10, ui.value);

  filter.frequency.value = cutoff;

  // setTimeout("drawCurve()", 50);
  drawCurve();
  var info = document.getElementById('cutoff-value');
  info.textContent = `cutoff = ${cutoff.toFixed(1)} Hz`;
}

function qHandler(event, ui) {
  var q = new Number(ui.value);
  filter.Q.value = q;
  // setTimeout("drawCurve()", 50);
  drawCurve();
  var info = document.getElementById('Q-value');
  info.textContent = `Q = ${q.toFixed(3)}`;
}

function qHandlerdB(event, ui) {
  var q = new Number(ui.value);
  filter.Q.value = q;
  // setTimeout("drawCurve()", 50);
  drawCurve();
  var info = document.getElementById('Q-value');
  info.textContent = `Q = ${q.toFixed(3)} dB`;
}



function gainHandler(event, ui) {
  console.log(filter);

  var gain = new Number(ui.value);
  filter.gain.value = gain;
  // setTimeout("drawCurve()", 100);
  drawCurve();
  var info = document.getElementById('gain-value');
  info.textContent = `gain = ${gain.toFixed(3)} dB`;
}

function outputGainHandler(event, ui) {
  let g = new Number(ui.value);
  // The output gain is in dB; convert to linear.
  outputGain.gain.value = Math.pow(10, g / 20);

  let element = document.getElementById('Volume-value');
  element.textContent = `Volume = ${g.toFixed(2)} dB`;
}

function setFilterType(filterType) {
  filter.type = filterType;

  // Need to adjust the Q slider range depending on the filter type.  This
  // is because the Q parameter for the lowpass and highpass filters is in
  // dB, but linear for the others.

  switch (filterType) {
    case 'lowpass':
    case 'highpass':
      // Q is in dB for these filters.  And the gain is not used
      document.getElementById('QSlider').disabled = false;
      configureSlider('Q', q, -50, 50, qHandlerdB);
      document.getElementById('gainSlider').disabled = true;
      break;
    case 'bandpass':
    case 'notch':
    case 'allpass':
      // Q is linear and gain is not used
      document.getElementById('QSlider').disabled = false;
      configureSlider('Q', Math.max(q, 0), 0, 100, qHandler);
      document.getElementById('gainSlider').disabled = true;
      break;
    case 'peaking':
      // Q is linear, and gain is used
      document.getElementById('QSlider').disabled = false;
      configureSlider('Q', q, 0, 100, qHandler);
      document.getElementById('gainSlider').disabled = false;
      break;
    case 'lowshelf':
    case 'highshelf':
      // Q is not used, but gain is.
      document.getElementById('QSlider').disabled = true;
      document.getElementById('gainSlider').disabled = false;
      break;
    default:
      console.log(`Unhandled filter type: ${filterType}`);
      break;
  }
  drawCurve();
}

function animateCurve() {
  drawCurve();
  requestAnimationFrame(animateCurve);
}

function init() {
  $('<div id=\'tooltip\'><\/div>')
      .css({
        position: 'absolute',
        display: 'none',
        border: '1px solid #fdd',
        padding: '2px',
        'background-color': '#fee',
        opacity: 0.8
      })
      .appendTo('body');

  AudioContext = AudioContext || webkitAudioContext;
  // Setup biquad and other audio stuff
  context = new AudioContext();
  nyquist = context.sampleRate / 2;
  filter = context.createBiquadFilter();
  filter.type = 'bandpass';         // Bandpass
  filter.frequency.value = cutoff;  // cutoff
  filter.Q.value = q;
  filter.gain.value = gain;
  outputGain = context.createGain();
  var period = 2;
  var startTime = context.currentTime;

  //      for (var i = 0; i < 100; ++i) {
  //          filter.frequency.exponentialRampToValueAtTime(200.0,
  //          startTime + i*period);
  //          filter.frequency.exponentialRampToValueAtTime(10000.0,
  //          startTime + i*period + 0.5*period);
  //      }


  filter.connect(outputGain).connect(context.destination);

  addSlider('cutoff');
  addSlider('Q');
  addSlider('gain');
  addSlider('Volume');
  // Default values for the sliders.  These may get reconfigured when the
  // selected filter type changes.

  // The cutoff slider is in log10 units ranging from lowestFrequency
  // to Nyquist.  Set set the limits to the log of these values, after
  // normalizing them.
  configureSlider(
      'cutoff', Math.log10(cutoff), Math.log10(lowestFrequency),
      Math.log10(nyquist), cutoffHandler);
  configureSlider('Q', q, 0, 100, qHandler);
  configureSlider('gain', gain, -40.0, 40.0, gainHandler);
  configureSlider('Volume', 0, -10, 10, outputGainHandler);

  // The default (checked) button below is bandpass, so set up sliders for
  // bandpass.
  setFilterType('bandpass');

  window.addEventListener('resize', adjustSliderPositions);

  // Give audio process some time initialize itself.
  drawCurve();
  // animateCurve();
}
