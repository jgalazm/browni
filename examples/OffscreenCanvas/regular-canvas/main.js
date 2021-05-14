const simulationParameters = {
  skip: 1
}

const workerProgram = (canvas) => {
  const ds = 15 / 60;
  const extent = {
    xmin: -121,
    xmax: 25,
    ymin: -60,
    ymax: 85,
  };

  const scenario = {
    ...extent,
    waveWidth: parseInt((extent.xmax - extent.xmin) / ds),
    waveHeight: parseInt((extent.ymax - extent.ymin) / ds),
    earthquake: {
      ce: -12, //centroid N coordinate, e
      cn: -35,
      Mw: 9.0
    },
  };

  const nami = new Nami(scenario, { canvas, loop: true }, {
    modelStepDidFinish: (model, thisController) => {
      if (model.discretization.stepNumber % simulationParameters.skip === 0) {
        return false;
      }
      return true;
    }
  });

}

const cycleSpeed = () => {
  simulationParameters.skip = (simulationParameters.skip + 10) % 110;
  document.querySelector("#speed-counter").innerText = `Skipping: ${simulationParameters.skip-1} frames`;
}

const mainProgram = () => {
  const canvas = document.querySelector("#simulation-canvas");
  canvas.width = canvas.offsetWidth;
  canvas.height = canvas.offsetHeight;

  workerProgram(canvas);
}


