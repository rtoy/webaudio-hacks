var sampleRate = 48000.0;

var renderLengthSeconds = 1;
var pulseLengthSeconds = 0.1;
var pulseLengthFrames = pulseLengthSeconds * sampleRate;

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

function checkConvolvedResult(trianglePulse) {
    return function(event) {
        var renderedBuffer = event.renderedBuffer;

        referenceData = trianglePulse.getChannelData(0);
        renderedData = renderedBuffer.getChannelData(0);
    
        var success = true;
    
        drawCurve("#conv", renderedData, 0, renderLengthSeconds * sampleRate);
        drawCurve("#expected", referenceData, 0, renderLengthSeconds * sampleRate);

        var n = renderedBuffer.length;

        var maxDelta = 0;
        var valueAtMaxDelta = 0;
        signalDiff = new Array();
        for (var i = 0; i < referenceData.length; ++i) {
            var diff = renderedData[i] - referenceData[i];
            signalDiff[i] = diff;
            var x = Math.abs(diff);
            if (x > maxDelta) {
                maxDelta = x;
                valueAtMaxDelta = referenceData[i];
                maxDeltaIndex = i;
            }
        }

        drawCurve("#diff", signalDiff, 0, renderLengthSeconds * sampleRate);

        // Make sure that portion after convolved portion is totally silent.
        var isFinalPortionSilent = true;
        for (var i = referenceData.length + 128; i < renderedBuffer.length; ++i) {
            if (renderedData[i] != 0) {
                isFinalPortionSilent = false;
                // alert(i + ": renderedData[i] = " + renderedData[i]);
            }
        }
        
        if (isFinalPortionSilent) {
            testPassed("Rendered signal after tail of convolution is silent.");
        } else {
            testFailed("Rendered signal after tail of convolution should be silent.");
        }
        
        var maxDeviationFraction = maxDelta / valueAtMaxDelta;

        console.log("n = " + n + " : maxDelta = " + maxDelta + " : valueAtMaxDelta = " + valueAtMaxDelta + " : maxDeltaIndex = " + maxDeltaIndex + " : maxDeviationFraction = " + maxDeviationFraction);

        // for (var i = maxDeltaIndex - 4; i < maxDeltaIndex + 4; ++i) {
        // for (var i = 0; i < referenceData.length; ++i) {
        //     var d1 = renderedData[i + 128];
        //     var d2 = referenceData[i];
        //     console.log(i + ": " + d1 + ", " + d2);
        // }
        
        if (success) {
            testPassed("Test signal was correctly delayed.");
        } else {
            testFailed("Test signal was not correctly delayed.");
        }

        finishJSTest();
    }
}
