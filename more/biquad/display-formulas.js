let filter;
function texifyNumber(number, options) {
    // Convert the number to a string.  If it is in scientific form,
    // replace with the appropriate teX version.
    if (options && !options.showUnity && number === 1) {
        if (options.addSign) {
            return Math.sign(number) < 0 ? "-" : "+";
        }
        return "";
    }
    if (options && options.addSign) {
        var s = Math.abs(number).toString();
        s = s.replace(/(e)(.*)/, "\\times 10^{$2}");
        s = (number < 0) ? "- " + s : "+" + s;
    } else {
        var s = number.toString();
        s = s.replace(/(e)(.*)/, "\\times 10^{$2}");
    }

    return s;
}

function digitalTermTeX(term) {
    var f = "\\frac{";
    if (term[0][1].length == 2) {
        // Linear term
        if (term[0][1][1] == 0) {
            f += "1";
        } else {
            f += "1 " + texifyNumber(term[0][1][1], {
                addSign: true
            }) + "z^{-1}";
        }
        f += "}{";
        f += "1 " + texifyNumber(term[1][1], {
            addSign: true
        }) + "z^{-1}";
        f += "}";
    } else {
        f += texifyNumber(term[0][0], {
            showUnity: false
        }) + "1";
        if (term[0][1][0] != 0) {
            f += texifyNumber(term[0][1][1], {
                addSign: true
            }) + "\\,z^{-1}";
        }
        if (term[0][1][2] != 0) {
            f += texifyNumber(term[0][1][2], {
                addSign: true
            }) + "\\,z^{-2}";
        }
        f += "}{";
        f += "1 " + texifyNumber(term[1][1], {
            addSign: true
        }) + "\\,z^{-1} ";
        f += texifyNumber(term[1][2], {
            addSign: true
        }) + "\\,z^{-2}";
        f += "}";
    }

    return f;
}

let filterType = "lowpass";
function setFilterType(type) {
    filterType = type;
}

function calc() {
    let sampleRate = document.getElementById("samplerate").value;
    let freq = document.getElementById("frequency").value;
    let Q = document.getElementById("Q").value;
    let gain = document.getElementById("gain").value;

    console.log("filterType " + filterType);
    console.log("sampleRate " + sampleRate);
    console.log("freq = " + freq);
    console.log("Q = " + Q);
    console.log("gain = " + gain);

    filter = createFilter(filterType, freq / sampleRate, Q, gain);
    console.log(filter);

    let term = [[1, [filter.b0, filter.b1, filter.b2]],
                 [1, filter.a1, filter.a2]];
    let formula = digitalTermTeX(term);
    console.log(formula);

    MathJax.typesetPromise().then(() => {
	    const eqn = document.querySelector("#eqn");
	    eqn.innerHTML = "$$" + formula + "$$";
	    MathJax.typesetPromise();
	}).catch((err) => console.log(err.message));
	    
}
