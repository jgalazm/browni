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
    loop:true,
    stopTime: 10,
    displayWidth: 512,
    displayHeight: 512,
};
var play, video1, video2, mediaRecorder;
var  recordedBlobs = [];

lifeCycle = {
    dataWasLoaded: (model) => {
        document.body.appendChild(model.canvas);
        model.canvas.style = 'position:absolute;left:-512px;';
        
        video1 = document.getElementById('videoTarget');
        var stream = model.canvas.captureStream();
        video1.srcObject = stream;
    
    }
}

setTimeout(function(){
    let thismodel = new NAMI.app(data, output, lifeCycle);
}, 2000);