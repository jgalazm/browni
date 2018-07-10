let Controller = function(model,data, lifeCycle){
    let thisController = {};

    
    let controllerSimulationDidFinish, modelStepDidFinish, iterationDidFinish;

    let paused = false;

    if( typeof lifeCycle !== 'undefined'){
        if( typeof lifeCycle.controllerSimulationDidFinish !== 'undefined'){
            controllerSimulationDidFinish = lifeCycle.controllerSimulationDidFinish;
        }
        else{
            controllerSimulationDidFinish = ()=>{};
        }

        if( typeof lifeCycle.modelStepDidFinish !== 'undefined'){
            modelStepDidFinish = lifeCycle.modelStepDidFinish;
        }
        else{
            modelStepDidFinish = ()=>{};
        }

        if( typeof lifeCycle.iterationDidFinish !== 'undefined'){
            iterationDidFinish = lifeCycle.iterationDidFinish;
        }
        else{
            iterationDidFinish = (model, controller, animate)=>{  requestAnimationFrame(animate); };
        }        
    }

    let downloadGridArray = function(array, time){
        
        array = array.toString();
        
        array = time.toString()+'\n'+array;

        var outputData = new Blob([array], { type: 'text/csv' }); 
        var outputURL = URL.createObjectURL(outputData);
        var link = document.createElement('a');
            
        link.href =  outputURL;
        
        link.download = 'tlab2D';

        link.click();  
    }

    let downloadPOIs = function(pois){
        var exportString = "#seconds";
        Object.keys(pois).forEach(function(poi){
            exportString +=','+poi+'(m)';
        });
        exportString+='\n';        

        
        for(let i=0;i< pois[Object.keys(pois)[0]].time.length; i++){
            exportString += pois[Object.keys(pois)[0]].time[i].toFixed(8).toString();
            
            
            Object.keys(pois).forEach(function(poi){
                exportString += ','+pois[poi].surface[i].toFixed(8).toString();
            });

            exportString += '\n';
        }
        
        exportString = 'data:text/csv;charset=utf-8,'+exportString;   
        exportString = encodeURI(exportString);     

        let link = document.createElement("a");
        link.href =  exportString;
        link.download = "tseries";
        link.click();   

        }

    let animate = () => {

        
        
        if(model.discretization.stepNumber == 0){
            if(lifeCycle.modelSimulationWillStart)
            lifeCycle.modelSimulationWillStart(model, thisController);
        }

        if(paused){
            requestAnimationFrame(animate);
            return;
        }

        /********************* */
        let exitLoop = true;
        while(exitLoop){
            model.runSimulationStep();                
                
            exitLoop = modelStepDidFinish(model, thisController);
        }

        model.displayPColor();
        /********************* */


        if(model.currentTime<data.stopTime){
            
            iterationDidFinish(model, thisController, animate);
            
        }
        else if(data.loop){

            model.newEarthquake = model.earthquake;
            iterationDidFinish(model, thisController, animate);

        }
        else{

            controllerSimulationDidFinish(model, thisController);
            
        }
    
    }

    thisController = {
        animate,
        togglePause: ()=>{
            paused = !paused;
        },
        set paused(value){
            paused = value;
        },
        get paused(){
            return paused;
        },
        downloadCurrentGridHeights: ()=> {
            downloadGridArray(model.currentGridHeights,
                model.discretization.stepNumber*model.discretization.dt);
        },
        downloadMaximumHeights: () =>{
            downloadGridArray(model.currentMaximumHeights,
                model.discretization.stepNumber*model.discretization.dt);
        },
        downloadArrivalTimes: () =>{
            downloadGridArray(model.currentArrivalTimes,
                model.discretization.stepNumber*model.discretization.dt);
        },
        downloadAllPois: () =>{
            downloadPOIs(model.pois);
        },
        get stopTime(){
            return data.stopTime;
        }
    };

    return thisController;
}

export {Controller};