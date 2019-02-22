import {Controller} from './Controller';
import {Model} from './Model';
import Reader from './Reader';

let Nami = function(data, output, lifeCycle){

    let model, controller;

    let init = () => {
        
        this.model = new Model(data, output);

        if (lifeCycle.dataWasLoaded !== undefined){
            lifeCycle.dataWasLoaded(this.model);
        }

        this.controller = new Controller(this.model, output, lifeCycle);

        this.controller.animate();
        
    }

    const reader = new Reader(data, init);

    

}

export default Nami;
