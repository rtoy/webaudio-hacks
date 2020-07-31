//
// Returns the filter response H(e^(jw)). coef is a dictionary of the
// coefficients of the digital filter.
//
//
function H(omegas, coef) {
  let {b0, b1, b2, a1, a2} = coef;

  let mag = [];
  let phase = [];
  omegas.forEach((omega) => {
    let w = Math.PI * omega;
    let s = Math.sin(w);
    let c = Math.cos(w);
    let s2 = Math.sin(2 * w);
    let c2 = Math.cos(2 * w);

    let topReal = b0 + b1 * c + b2 * c2;
    let topImag = -(b1 * s + b2 * s2);
    let botReal = 1 + a1 * c + a2 * c2;
    let botImag = -(a1 * s + a2 * s2);

    if (1) {
      // Perform complex division
      let denom = Math.pow(Math.hypot(botReal, botImag), 2);

      let responseReal = (topReal * botReal + topImag * botImag) / denom;
      let responseImag = (topImag * botReal - topReal * botImag) / denom;

      mag.push(
          [omega, 20 * Math.log10(Math.hypot(responseReal, responseImag))]);
      phase.push([omega, Math.atan2(responseImag, responseReal)]);
    } else {
      let responseMag =
          Math.hypot(topReal, topImag) / Math.hypot(botReal, botImag);
      let responsePhase =
          Math.atan2(topImag, topReal) - Math.atan2(botImag, botReal);
      mag.push([omega, 20 * Math.log10(responseMag)]);
      phase.push([omega, responsePhase]);
    }
  });
  return [mag, phase];
}

function getResponse(filter, sampleRate) {
  const nyquist = sampleRate / 2;

  // How many samples to take for the frequency axis.
  const steps = 1000;

  // Array of frequency samples.  The values are normalized
  // frequencies from 0 to 1 where 1 represents Nyquist.
  let omega = new Array(steps);

  // Lowest frequency we want to show.  Should be a power of 10 so the
  // graph starts at a nice power of 10 location.
  let lowestFrequency = document.getElementById('lowest-freq').value;;

  // Logarithmically sample between the lowest frequency and Nyquist.
  let loFreq = Math.log10(lowestFrequency);
  let hiFreq = Math.log10(nyquist);
  let delta = (hiFreq - loFreq) / steps

  for (let k = 0; k < steps; ++k) {
      let f = loFreq + k * delta;
      omega[k] = Math.pow(10, f) / nyquist;
  }

  let response = H(omega, filter);

  return response;
}

let logaxis = true;

function setLogAxis(isLog) {
    logaxis = isLog;
    calc();
}

function plotResponse(filterType, filter, sampleRate) {
  const nyquist = sampleRate / 2;
  let response = getResponse(filter, sampleRate);
  let magResponse = response[0];
  let phaseResponse = response[1];

  // Convert the normalized frequencies back to Hz.  And convert radians to
  // degrees
  const radToDeg = 180 / Math.PI;
  for (let k = 0; k < magResponse.length; ++k) {
    magResponse[k][0] = magResponse[k][0] * nyquist;
    phaseResponse[k][0] = phaseResponse[k][0] * nyquist;
    phaseResponse[k][1] = radToDeg * phaseResponse[k][1];
  }

  let ymin;
  let ymax;
  if (filterType == 'allpass') {
    ymin = -1;
    ymax = 1;
  }
  
      
  //$.plot("#graph"), [response[0], response[1]]);
  let graph = $('#graph');

  let plot = $.plot(
      graph,
      [
        {
          data: magResponse,
          label: 'Mag (dB)',
          lines: {lineWidth: 2},
          color: 'red'
        },
        {
          data: phaseResponse,
          label: 'Phase (deg)',
          lines: {lineWidth: 2},
          color: 'green',
          yaxis: 2
        }
      ],
      {
        grid: { hoverable: true},
        legend: {
          position: 'ne',
          show: true,
          // container: document.getElementById("legendContainer")
        },
        xaxis: {
          mode: logaxis ? 'log' : null,
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
              tickDecimals: 1
            // ticks : tickScale
          },
          {
            axisLabel: 'Phase (deg)',
            // align if we are to the right
            alignTicksWithAxis: 1,
            position: 'right',
            // tickFormatter: degFormatter,
            min: -180,
            max: 180,
            autoScale: 'none',
            showMinorTicks: true,
            //  ticks: [-180, -135, -90, -45, 0, 45, 90, 135, 180]
          }
        ],
      });
  graph.resize(() => {
    $('.message')
        .text(
            'Placeholder is now ' + $(this).width() + 'x' + $(this).height() +
            ' pixels');
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
	      $('#tooltip').html(
		  `${item.series.label} at ${x} = ${y}`)
		  .css({top: item.pageY+5, left: item.pageX+5})
		  .fadeIn(200);
	  } else {
	      $('#tooltip').stop().hide();
	  }
      });
  $('#graph').bind('plothovercleanup', (event, pos, item) => {
	  $('#tooltip').hide();
      });
}
