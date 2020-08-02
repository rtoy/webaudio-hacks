function addSlider(name) {
  let controls = document.getElementById('controls');

  let divName = name + 'Slider';

  let sliderText = `<div class="row">
       <div class="column">
         <input class="slider-bar "id="${divName}" type="range"
          min="0" max="0" step="0.01" value="0">
       </div>
       <div class="column">
         <div class="slider-text" id="${name}-value">
           ${name}
         </div.
       </div>
     </div>
     `;

  controls.innerHTML = controls.innerHTML + sliderText;
}

function configureSlider(name, value, min, max, handler) {
  let divName = name + 'Slider';

  let slider = document.getElementById(divName);

  slider.min = min;
  slider.max = max;
  slider.value = value;

  // Use oninput so the biquad changes as we move the sliders around
  // instead of waiting until we stop moving the slider.  To go back
  // to the old behavior, use onchange instead of oninput.
  slider.oninput = function() {
    handler(0, this);
  };
}
