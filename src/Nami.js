import { Controller } from "./Controller";
import { Model } from "./Model/Model";
import Reader from "./Reader/Reader";

let Nami = function(data, output, lifeCycle) {
  let model, controller;
  let dataWasParsed = () =>{};
  if(lifeCycle && lifeCycle.dataWasParsed) dataWasParsed = lifeCycle.dataWasParsed;

  let dataWasLoaded = (model) =>{
    document.body.appendChild(model.canvas);
  };
  if(lifeCycle && lifeCycle.dataWasLoaded) dataWasLoaded = lifeCycle.dataWasLoaded;

  let init = newData => {
    dataWasParsed(newData);

    this.model = new Model(newData);

    dataWasLoaded(this.model);

    this.controller = new Controller(this.model, newData, lifeCycle);

    this.controller.animate();
  };

  const newData = new Reader(data, output);
  Promise.all([newData.bathymetry.array, newData.initialCondition, newData.slab]).then(
    values => {
      const [bathyArray, earthquake, slab] = values;

      newData.bathymetry.array = bathyArray;

      newData.earthquake = earthquake;

      newData.slab = slab;
      
      init(newData);
    }
  );
};

export default Nami;
