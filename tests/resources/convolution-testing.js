var sampleRate = 48000.0;

var referenceData;
var renderedData;
var signalDiff;

function createSquarePulseBuffer(context, sampleFrameLength) {
    var audioBuffer = context.createBuffer(1, sampleFrameLength, context.sampleRate);

    var n = audioBuffer.length;
    var data = audioBuffer.getChannelData(0);

    for (var i = 0; i < n; ++i)
        data[i] = 1;

    return audioBuffer;
}

// The triangle buffer holds the expected result of the convolution.
// It linearly ramps up from 0 to its maximum value (at the center) then linearly ramps down to 0.
// The center value corresponds to the point where the two square pulses overlap the most.
function createTrianglePulseBuffer(context, sampleFrameLength) {
    var audioBuffer = context.createBuffer(1, sampleFrameLength, context.sampleRate);

    var n = audioBuffer.length;
    var halfLength = n / 2;
    var data = audioBuffer.getChannelData(0);
    
    var maxValue = halfLength;

    for (var i = 0; i < halfLength; ++i)
        data[i] = i + 1;

    for (var i = halfLength; i < n; ++i)
        data[i] = n - i - 1;

    return audioBuffer;
}

function log10(x) {
  return Math.log(x)/Math.LN10;
}

function linearToDecibel(x) {
  return 20*log10(x);
}

// Verify that the rendered result is very close to the reference
// triangular pulse.
function checkTriangularPulse(rendered, reference) {
    var match = true;
    var maxDelta = 0;
    var valueAtMaxDelta = 0;
    var maxDeltaIndex = 0;

    for (var i = 0; i < reference.length; ++i) {
        var diff = rendered[i] - reference[i];
        var x = Math.abs(diff);
        if (x > maxDelta) {
            maxDelta = x;
            valueAtMaxDelta = reference[i];
            maxDeltaIndex = i;
        }
    }

    // allowedDeviationFraction was determined experimentally.  It
    // is the threshold of the relative error at the maximum
    // difference between the true triangular pulse and the
    // rendered pulse.
    var allowedDeviationDecibels = -132.4;
    var maxDeviationDecibels = linearToDecibel(maxDelta / valueAtMaxDelta);

    if (maxDeviationDecibels <= allowedDeviationDecibels) {
        testPassed("Triangular portion of convolution is correct.  Max deviation = " + maxDeviationDecibels + " dB at " + maxDeltaIndex);
    } else {
        testFailed("Triangular portion of convolution is not correct.  Max deviation = " + maxDeviationDecibels + " dB at " + maxDeltaIndex);
        match = false;
    }

    return match;
}        

// Verify that the rendered data is close to zero for the first part
// of the tail.
function checkTail1(data, reference, breakpoint) {
    var isZero = true;
    var tail1Max = 0;

    for (var i = reference.length; i < reference.length + breakpoint; ++i) {
        var mag = Math.abs(data[i]);
        if (mag > tail1Max) {
            tail1Max = mag;
        }
    }

    // Let's find the peak of the reference (even though we know a
    // priori what it is).
    var refMax = 0;
    for (var i = 0; i < reference.length; ++i) {
      refMax = Math.max(refMax, Math.abs(reference[i]));
    }

    // This threshold is experimentally determined by examining the
    // value of tail1MaxDecibels.
    var threshold1 = -129.7;

    var tail1MaxDecibels = linearToDecibel(tail1Max/refMax);
    if (tail1MaxDecibels <= threshold1) {
        testPassed("First part of tail of convolution is sufficiently small: " + tail1MaxDecibels + " dB");
    } else {
        testFailed("First part of tail of convolution is not sufficiently small: " + tail1MaxDecibels + " dB");
        isZero = false;
    }

    return isZero;
}

// Verify that the second part of the tail of the convolution is
// exactly zero.
function checkTail2(data, reference, breakpoint) {
    var isZero = true;
    var tail2Max = 0;
    // For the second part of the tail, the maximum value should be
    // exactly zero.
    var threshold2 = 0;
    for (var i = reference.length + breakpoint; i < data.length; ++i) {
        if (Math.abs(data[i]) > 0) {
            isZero = false; 
            break;
        }
    }

    if (isZero) {
        testPassed("Rendered signal after tail of convolution is silent.");
    } else {
        testFailed("Rendered signal after tail of convolution should be silent.");
    }

    return isZero;
}

function drawDiff(actual, reference, start, stop) {
    signalDiff = new Array();
    for (var k = 0; k < reference.length; ++k) {
        signalDiff[k] = actual[k] - reference[k];
    }
    for (var k = reference.length; k < actual.length; ++k) {
        signalDiff[k] = actual[k];
    }

    drawCurve("#diff", signalDiff, start, stop);
}

function checkConvolvedResult(trianglePulse) {
    return function(event) {
        var renderedBuffer = event.renderedBuffer;

        referenceData = trianglePulse.getChannelData(0);
        renderedData = renderedBuffer.getChannelData(0);
    
        var success = true;
    
        drawCurve("#conv", renderedData, 0, renderLengthSeconds * sampleRate);
        drawCurve("#expected", referenceData, 0, renderLengthSeconds * sampleRate);

        drawDiff(renderedData, referenceData, 0, renderLengthSeconds * sampleRate);
        
        // Verify the triangular pulse is actually triangular.

        success = success && checkTriangularPulse(renderedData, referenceData);
        
        // Make sure that portion after convolved portion is totally
        // silent.  But round-off prevents this from being completely
        // true.  At the end of the triangle, it should be close to
        // zero.  If we go farther out, it should be even closer and
        // eventually zero.

        // For the tail of the convolution (where the result would be
        // theoretically zero), we partition the tail into two
        // parts.  The first is the at the beginning of the tail,
        // where we tolerate a small but non-zero value.  The second part is
        // farther along the tail where the result should be zero.
        
        // breakpoint is the point dividing the first two tail parts
        // we're looking at.  Experimentally determined for the given
        // pulse length.
        var breakpoint = 9800;

        success = success && checkTail1(renderedData, referenceData, breakpoint);
        
        success = success && checkTail2(renderedData, referenceData, breakpoint);
        
        if (success) {
            testPassed("Test signal was correctly convolved.");
        } else {
            testFailed("Test signal was not correctly convolved.");
        }

        finishJSTest();
    }
}
