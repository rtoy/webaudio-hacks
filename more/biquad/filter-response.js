//
// Returns the filter response H(e^(jw)). coef is a dictionary of the
// coefficients of the digital filter.
//
// 
function H(omegas, coef) {
  let {b0,b1,b2,a1,a2} = coef;

  let mag = [];
  let phase = [];
  omegas.forEach((omega) => {
      let s = Math.sin(omega);
      let c = Math.cos(omega);
      let s2 = Math.sin(2*omega);
      let c2 = Math.cos(2*omega);

      let topReal = b0 + b2*c + b2*c2;
      let topImag = -(b1*s + b2*s2);
      let botReal = 1 + a1*c + a2 * c2;
      let botImag = -(a1*s + a2*s2);

      if (0) {
	// Perform complex division
	let denom = Math.pow(Math.hypot(botReal, botImag), 2);

	let responseReal = (topReal * botReal + topImag * botImag) / denom;
	let responseImag = (topImag * botReal - topReal * botImag) / denom;

	mag.push([omega, 20*Math.log10(Math.hypot(responseReal, responseImag))]);
	phase.push([omega, Math.atan2(responseImag, responseReal)]);
      } else {
	let responseMag = Math.hypot(topReal, topImag) / Math.hypot(botReal, botImag);
	let responsePhase = Math.atan2(topImag, topReal) - Math.atan2(botImag, botReal);
	mag.push([omega, 20*Math.log10(responseMag)]);
	phase.push([omega, responsePhase]);
      }
    });
  return [mag, phase];
}

function getResponse(filter, sampleRate) {
  const noctaves = 8;
  const nyquist = sampleRate / 2;
  // Just uniformly sample from 0 to pi.
  const steps = 1000;
  let omega = new Array(steps);

  for (let k = 0; k < steps; ++k) {
    let f = k / steps;
    // Conver to log frequency scale (octaves)
    f = Math.pow(2, noctaves * (f - 1));  
    omega[k] = f;
  }

  let response = H(omega, filter);

  return response;
}

function plotResponse(filter, sampleRate) {
  const nyquist = sampleRate / 2;
  let response = getResponse(filter, sampleRate);
  let magResponse = response[0];
  let phaseResponse = response[1];

  // Convert the normalized frequencies back to Hz.  And convert radians to degrees
  const radToDeg = 180 / Math.PI;
  for (let k = 0; k < magResponse.length; ++k) {
    magResponse[k][0] = magResponse[k][0] * nyquist;
    phaseResponse[k][0] = phaseResponse[k][0] * nyquist;
    phaseResponse[k][1] = radToDeg * phaseResponse[k][1];
  }
  //$.plot("#graph"), [response[0], response[1]]);
  $.plot($("#graph"),
         [{
	     data : magResponse,
	     label: "Mag (dB)",
	     lines: {linewidth: 3}
	   },
          {
	     data : phaseResponse,
	     label: "Phase (deg)",
	     lines: {linewidth: 3},
	     yaxis: 2
	   }],
         {
	   xaxis: { mode: "log", ticks: 10, showMinorTicks: true },
	     //xaxes: [ { ticks : tickScale } ],
	     yaxes: [ { //tickFormatter: dBFormatter,
	       //min: magmin,
	       //max: magmax,
	       //ticks : tickScale
	     },
	       {
		 // align if we are to the right
		 alignTicksWithAxis: position = "right" ? 1 :
		 null,
		   position: position,
		   //tickFormatter: degFormatter,
		   //min: phasemin,
		   //max: phasemax,
		   //ticks : tickScale
		   }],
             legend: { position: 'ne' }
         }
    );
}

