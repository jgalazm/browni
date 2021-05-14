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
  simulationParameters.skip = simulationParameters.skip == 1 ? 51 : 1;
  document.querySelector("#speed-counter").innerText = `Skipping: ${simulationParameters.skip-1} frames`;
}


let timeEl;
let t1 = Date.now();
let t2;
function tickTimer() {
  t2 = Date.now();
  const diff = t2 - t1;
  t1 = t2;
  timeEl.innerText = `approx ${(1000/diff).toFixed(1)} fps`;
  window.timerRafId = window.requestAnimationFrame(tickTimer);
}

const mainProgram = () => {
  const canvas = document.querySelector("#simulation-canvas");
  canvas.width = canvas.offsetWidth;
  canvas.height = canvas.offsetHeight;

  workerProgram(canvas);
  timeEl = document.getElementById('timer');
  requestAnimationFrame(tickTimer);
}


