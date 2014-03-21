// Support for both prefixed and unprefixed AudioContext and OfflineAudioContext.
// We don't support other legacy names.
if (window.hasOwnProperty('webkitAudioContext') &&
    !window.hasOwnProperty('AudioContext')) {
    window.AudioContext = webkitAudioContext;
    window.OfflineAudioContext = webkitOfflineAudioContext;
}
