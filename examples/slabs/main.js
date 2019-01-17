let w = 501;
let h = 501;

let slab = {
    'x': [-7.5, 7.5],
    'y': [0,0],
    'depth': [22,-40],
    'dip': [13, 13],
    'strike': [-45, 45]
}

let data = {
    xmin : -15,
    xmax :  15,
    ymin :  -15,
    ymax : 15,
    waveWidth: w,
    waveHeight: h,
    coordinates: 'spherical',
    bathymetry: [[1000,1000],[1000,1000]],
    earthquake: [{
        depth: 22900,
        strike: 17,
        dip: 13.0,
        rake: 108.0,
        U3: 0.0,
        cn: 0,   //centroid N coordinate, e
        ce: -5,
        Mw: 9.5,
        reference: 'mid bottom'
    }],
    slab: slab
}

let output = {
    displayWidth:  w,
    displayHeight: h,
    stopTime: 60*60*15,
    displayOption: 'heights',

};

let expectedOutput = 300;
let lifeCycle = {
    dataWasLoaded: (model)=>{
        document.body.appendChild(model.canvas);
    },
    

    modelStepDidFinish: (model, controller) =>{
        if(model.discretization.stepNumber % 10 === 0){
            console.log(model.discretization.stepNumber);
        }

        if(model.discretization.stepNumber>800){
            const ce = -model.earthquake[0].ce;
            model.earthquake = [{
                depth: 22900,
                rake: 108.0,
                U3: 0.0,
                cn: 0,   //centroid N coordinate, e
                ce: ce,
                Mw: 9.5,
                reference: 'mid bottom'
            }]
        }
        
        if(model.discretization.stepNumber % 10 !== 0){
            return true;
        }
        return false;
    },

}

let thisApp = new NAMI.app(data, output, lifeCycle);
