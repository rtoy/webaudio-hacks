let src = null;
let f = null;
let audioOn = false;

function createBasicGraph(sampleRate) {
  context = new AudioContext({sampleRate: sampleRate});
  mainGain = new GainNode(context, {gain: .5});
  mainGain.connect(context.destination);
}

async function enableAudio() {
  if (audioOn) {
    audioOn = false;
    await context.suspend();
  } else {
    audioOn = true;
    await context.resume();
    playAudio(filterType, filterCoef);
  }
}

function playAudio(filterType, filterCoef) {
  if (!audioOn) {
    return;
  }

  let sampleRate = document.getElementById('samplerate').value;
  if (sampleRate != context.sampleRate) {
    // Need to rebuild the graph
    createBasicGraph(sampleRate);
  }

  if (filterType == 'lowshelfq' || filterType == 'highshelfq') {
    // WebAudio currently doesn't have these types of filters, so we
    // fake it with an IIRFilter.  We don't need automation, so this
    // is fine.
    f = new IIRFilterNode(context, {
	feedforward: [filterCoef.b0, filterCoef.b1, filterCoef.b2],
	  feedback: [1, filterCoef.a1, filterCoef.a2]
	  });
  } else {
    // Set up regular BiquadFilter
    let freq = document.getElementById('frequency').value;
    let Q = document.getElementById('Q').value;
    let gain = document.getElementById('gain').value;

    f = new BiquadFilterNode(context, {
	type: filterType,
	  frequency: freq,
	  Q: Q,
	  gain: gain
	  });
  }

  f.connect(mainGain);

  if (src) {
    src.stop();
    src = null;
  }
  src = new OscillatorNode(context, {type: 'sawtooth', frequency: 880});
  src.connect(f);
  src.start();
}
