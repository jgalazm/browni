import { Controller } from "./Controller";
import { Model } from "./Model/Model";
import Reader from "./Reader/Reader";

let Nami = function(data, output, lifeCycle) {
  let model, controller;

  let init = newData => {
    if (lifeCycle.dataWasParsed) {
      lifeCycle.dataWasParsed(newData);
    }

    this.model = new Model(newData, output);

    if (lifeCycle.dataWasLoaded !== undefined) {
      lifeCycle.dataWasLoaded(this.model);
    }

    this.controller = new Controller(this.model, output, lifeCycle);

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
