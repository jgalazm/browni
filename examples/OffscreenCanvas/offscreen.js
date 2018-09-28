importScripts('../../build/nami.js')

onmessage = function(evt) {
    let {canvas} = evt.data;
    
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
        ymax: 10,
        canvas: canvas
    }
    
    let output = {
        stopTime: 10,
        displayWidth: 512,
        displayHeight: 512,
    };
    
    let recordedBlobs = [];
    
    lifeCycle = {
        dataWasLoaded: (model) => {
            
            console.log(canvas);
            
            // offscreen.convertToBlob().then(function(blob) {
            //     let video = document.getElementById('video');
            //     var superBuffer = new Blob(recordedBlobs, {type: 'video/webm'});
            //     video.src = window.URL.createObjectURL(superBuffer);
            //     console.log(blob);
            // });
    
        },
    
        modelStepDidFinish: (model, controller) => {
            console.log(model.discretization.stepNumber);
            return false;
        }
    }
    
    let thismodel = new NAMI.app(data, output, lifeCycle);
};

