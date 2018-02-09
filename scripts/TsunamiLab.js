function TsunamiLab(inputParameters, outputParameters) {

    let bathymetry = inputParameters.bathymetry; // this could be a setter
    bathymetry.texture = new THREE.Texture(bathymetry.img);
    bathymetry.texture.needsUpdate = true;

    let finiteFaultParameters = inputParameters.finiteFaultParameters;

    let numericalMeshResolution = inputParameters.numericalMeshResolution;
    
    let CFL = inputParameters.CFL;

    let coriolisFlag = inputParameters.coriolisFlag;

    let pcolorMapResolution = outputParameters.pcolorMapResolution;
    
    let colormap = outputParameters.colormap;
    
    let pois = outputParameters.pois;

    let poisTimeStep = outputParameters.poisTimeStep;

    let gridsTimeStep = outputParameters.gridsTimeStep;

    let simulationStepsBetweenFrames = outputParameters.simulationStepsBetweenFrames;

    
    // initialize model 

    let modelParameters = {
        bathymetry,
        finiteFaultParameters,
        numericalMeshResolution,
        CFL,
        coriolisFlag,
        pcolorMapResolution,
        colormap,
        pois
    }

    let model = this.Model(modelParameters);



    let viewParameters = {
        canvas: model.getCanvas()
    };

    let view = this.View(viewParameters);

    let controllerParameters = {
        model,
        view,
        poisTimeStep,
        gridsTimeStep,
        simulationStepsBetweenFrames,
    }

    let controller = this.Controller(controllerParameters);

    return{
        controller,
        model,
        view
    }

}

