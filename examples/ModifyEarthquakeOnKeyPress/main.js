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
    waveWidth: 512,
    waveHeight: 512,
    xmin: -10,
    xmax: 10,
    ymin: -10,
    ymax: 10
}

let output = {
    stopTime: 1,
    displayWidth: 512,
    displayHeight: 512,
};

let recordedBlobs = [];

lifeCycle = {
    dataWasLoaded: (model) => {
        
        document.body.appendChild(model.canvas);
    },

    modelStepDidFinish: (model, controller) => {
        console.log(model.discretization.stepNumber, model.currentTime);
        return false;
    }
}

let thisApp = new NAMI.app(data, output, lifeCycle);


window.onkeydown = (ev)=>{
    if(ev.key === 'ArrowUp' || ev.key === 'ArrowDown'  ){

        const sign = ev.key === 'ArrowUp'? +1: -1;
        
        const currentDepth = thisApp.model.earthquake[0].depth;    
        const newDepth = Math.max(currentDepth + sign,1);
        const newSubfault = Object.assign(thisApp.model.earthquake[0], {depth: newDepth});
        const newEarthquake = [newSubfault];
        thisApp.model.earthquake = newEarthquake;

        const depthElement =  document.getElementById('depth');
        depthElement.innerHTML = `Current depth: ${newDepth} (m)`;
        console.log(depthElement.textContext);
        
        return;
    }
    
    console.log(ev.key);
};