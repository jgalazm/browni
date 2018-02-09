TsunamiLab.prototype.Controller = function(controllerParameters){
    
    // user parameters

    let model = controllerParameters.model,

        view = controllerParameters.view,
        
        poisTimeStep = controllerParameters.poisTimeStep,

        gridsTimeStep = controllerParameters.gridsTimeStep,
        
        simulationStepsBetweenFrames = controllerParameters.simulationStepsBetweenFrames;

    

    let runSimulationTimeStep = function(){
        
        
        let lastPOIsTime= model.poisTime[model.poisTime.length-1];

        let storePOIS = false;


        let lastGridsTime = model.gridsTime[model.gridsTime.length-1];

        let storeGrids = false;


        let nextIterationTime = model.simulationData.nextIterationTime;

        let currentIterationTime = model.simulationData.currentIterationTime;
        
        
        let oldCFL  = model.simulationData.CFL;



        if(lastPOIsTime + poisTimeStep < nextIterationTime){
            
            model.simulationData.timeStep = (lastPOIsTime + poisTimeStep - currentIterationTime)*60;

            storePOIS = true;
        
        }

        if(lastGridsTime + gridsTimeStep < nextIterationTime){
            
            model.simulationData.timeStep = (lastGridsTime + gridsTimeStep - currentIterationTime)*60;

            storeGrids = true;

            // if pois/grids output times are the same up to 0.01 seconds
            // then only modify time step once
            
            if(Math.abs(poisTimeStep - gridsTimeStep)>0.01/60){

                storePOIS = false;
            }
        }

        model.renderSimulation();

        // TODO: should add an option for rendering the pcolormap 
        // after skiping some iterations


        if(storePOIS){
            
            model.storePOISValues();
            
            model.simulationData.CFL = oldCFL;
        }

        if(storeGrids){
            
            model.store2DHeightValues();
            
            model.simulationData.CFL = oldCFL;

        }

    }

    let runAnimation = function(timeEnd,options){
        /*
            Runs the simulation and displays an animation
            until timeEnd (minutes) pass.
        */

        // TODO: hacer que termine justo en timeEnd

        let loop = function(){
            
            
            for(var i = 0; i < simulationStepsBetweenFrames; i++ ){
            
                if( model.simulationData.currentIterationTime <= timeEnd ){
                    
                    runSimulationTimeStep();

                    model.renderScreenVoid();
                    

                }
            }

            if( model.simulationData.currentIterationTime <= timeEnd ){

                requestAnimationFrame(loop);

            }
        };

        loop();
    }

    let getPOIs = function(){
        return {
            time: model.poisTime,
            pois: model.pois
        }
    }

    let getGrids = function(){
        return {
            time: model.gridsTime,
            surface: model.gridValues
        }
    }

    let getASCIIGrid = function(gridHeights){
        
        return model.convertGridToString(gridHeights);
    
    }
    
    let getASCIIPois = function(){
        var pois = getPOIs();
    
        // save pois to file

        var exportString = 'minutes';   
        
        Object.keys(pois.pois).forEach(function(poi){
            exportString +=','+poi;
        });

        exportString+='\n';
        
        for(var i=0;i< pois.time.length;i++){
        
            exportString += pois.time[i].toFixed(8).toString();


            Object.keys(pois.pois).forEach(function(poi){
                exportString += ','+pois.pois[poi].surface[i].toFixed(8).toString();
            });

            exportString += '\n';
        }
        exportString = 'data:text/csv;charset=utf-8,'+exportString;

        exportString = encodeURI(exportString);
        
        return exportString;
    }


    
    $('#poisdown').click(function(){
        
        let exportString = getASCIIPois();

        var link = document.createElement('a');

        link.href =  exportString;
        link.download = 'tseries';
        link.click();
        $(link).remove();

    });
    

    $('#gridDownload').click(function(){
        var grids = tsunamiLab.controller.getGrids();
        
        grids.surface.forEach(function(gridHeights,index){

            var outputURL = tsunamiLab.controller.getASCIIGrid(gridHeights);

            
            var link = document.createElement('a');
            
            link.href =  outputURL;
            
            link.download = 'tlab2D.'+(index).toString().padStart(5,'0');

            link.click();  
            $(link).remove();
        });
    });
    
    return {
        runSimulationTimeStep,
        runAnimation,
        getPOIs,
        getGrids,
        getASCIIGrid,
        getASCIIPois
    }
}