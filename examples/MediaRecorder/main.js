let data = {
    bathymetry: [[1,1,1],[1,1,1],[1,1,1]]
    ,
    earthquake: [{
        L: 5,
        W: 3,
        depth: 2,
        slip: 1.0,
        strike: 30.0,
        dip: 70.0,
        rake: -45.0,
        U3: 1.0,
        cn: 0,
        ce: 0,
        reference: 'center'
    }],
    coordinates: 'cartesian',
    waveWidth: 100,
    waveHeight: 100,
    xmin: -10,
    xmax: 10,
    ymin: -10,
    ymax: 10
}

let output = {
    stopTime: 10,
    displayWidth: 512,
    displayHeight: 512,
};
var play, video1, video2, mediaRecorder;
var  recordedBlobs = [];

lifeCycle = {
    dataWasLoaded: (model) => {
        document.body.appendChild(model.canvas);
        
        video1 = document.getElementById('videoSource');
        var stream = model.canvas.captureStream();
        video1.srcObject = stream;
        
        var interval = setInterval(function(){
            promise = video1.play();
            if (promise !== undefined) {
                promise.then(_ => {
                    clearInterval(interval);
                }).catch(error => {
                });
            }
        }, 100);
        
        
        var options = { mimeType: 'video/webm' };
        mediaRecorder = new MediaRecorder(stream, options);
        mediaRecorder.ondataavailable = handleDataAvailable;
        mediaRecorder.start(1000);
    
    
        function handleDataAvailable(event) {
          if (event.data && event.data.size > 0) {
            recordedBlobs.push(event.data);
          }
          video2 = document.getElementById('videoRecorded');
          var superBuffer = new Blob(recordedBlobs, {type: 'video/webm'});
          video2.src = window.URL.createObjectURL(superBuffer);
          video2.controls = true;
        }
    }
}

setTimeout(function(){
    let thismodel = new NAMI.app(data, output, lifeCycle);
}, 2000);