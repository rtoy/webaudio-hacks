// Support for both prefixed and unprefixed AudioContext and OfflineAudioContext.
// We don't support other legacy names.
if (!((typeof AudioContext === "function")
      || (typeof AudioContext === "object")
      || (typeof webkitAudioContext === "function")
      || (typeof webkitAudioContext === "object")
      )) {
  alert("Sorry! Web Audio not supported by this browser");
}

if (window.hasOwnProperty('webkitAudioContext') &&
    !window.hasOwnProperty('AudioContext')) {
    window.AudioContext = webkitAudioContext;
    window.OfflineAudioContext = webkitOfflineAudioContext;
}
