// Support for both prefixed and unprefixed AudioContext and OfflineAudioContext.
// We don't support other legacy names.
if (!((typeof webkitAudioContext === "function")
      || (typeof AudioContext === "function")
      || (typeof webkitAudioContext === "object")
      || (typeof AudioContext === "object"))) {
  alert("Sorry! Web Audio not supported by this browser");
}
if (window.hasOwnProperty('webkitAudioContext') &&
    !window.hasOwnProperty('AudioContext')) {
    window.AudioContext = webkitAudioContext;
    window.OfflineAudioContext = webkitOfflineAudioContext;
}
