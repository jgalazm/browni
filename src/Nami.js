import {Controller} from './Controller';
import {Model} from './Model';
import Reader from './Reader/Reader';

let Nami = function(data, output, lifeCycle){

    let model, controller;

    let init = () => {
        
        if(lifeCycle.dataWasParsed){
            lifeCycle.dataWasParsed(data);
        }
        
        this.model = new Model(data, output);

        if (lifeCycle.dataWasLoaded !== undefined){
            lifeCycle.dataWasLoaded(this.model);
        }

        this.controller = new Controller(this.model, output, lifeCycle);

        this.controller.animate();
        
    }

    const newData = new Reader(data);
    Promise.all([newData.bathymetry, newData.initialCondition]).then( values=>{

        const [bathyArray, earthquake] = values;

        data.bathymetry = {
            array : bathyArray 
        }

        data.earthquake = earthquake;
        init();
    })

    

}

export default Nami;
