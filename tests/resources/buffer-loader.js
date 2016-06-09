function BufferLoader(context, urlList, callback) {
    this.context = context;
    this.urlList = urlList;
    this.onload = callback;
    this.bufferList = new Array();
    this.loadCount = 0;
}

BufferLoader.prototype.loadBuffer = function(url, index) {
    // Load buffer asynchronously
    var request = new XMLHttpRequest();
    request.open("GET", url, true);
    request.responseType = "arraybuffer";

    var loader = this;

    request.onload = function() {
        loader.context.decodeAudioData(
            request.response,
            function (decodedAudio) {
                try {
                    loader.bufferList[index] = decodedAudio;
                    if (++loader.loadCount == loader.urlList.length)
                        loader.onload(loader.bufferList);
                } catch(e) {
                    alert('BufferLoader: unable to load buffer' + index);
                }
            },
            function () {
                alert('error decoding file data: ' + url);
            });
    }

    request.onerror = function() {
        alert('BufferLoader: XHR error');        
    }

    request.send();
}

BufferLoader.prototype.load = function() {
    for (var i = 0; i < this.urlList.length; ++i)
        this.loadBuffer(this.urlList[i], i);
}
