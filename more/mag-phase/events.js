function addSlider(name) {
  let controls = document.getElementById('controls');

  let divName = name + 'Slider';


  let sliderText = '<div style="width:500px; height:20px;"> <input id="' +
      divName + '" ' +
      'type="range" min="0" max="1" step="0.01" value="0" style="height: 20px; width: 450px;"> <div id="' +
      name + '-value" style="position:relative; left:30em; top:-18px;">' +
      name + '</div> </div> <br>  ';

  controls.innerHTML = controls.innerHTML + sliderText;
}

function configureSlider(name, value, min, max, handler) {
  // let controls = document.getElementById("controls");
  //

  let divName = name + 'Slider';

  let slider = document.getElementById(divName);

  slider.min = min;
  slider.max = max;
  slider.value = value;
  slider.onchange = function() {
    handler(0, this);
  };
}

