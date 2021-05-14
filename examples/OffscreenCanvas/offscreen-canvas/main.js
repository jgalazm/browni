
const workerProgram = () => {
  importScripts("http://localhost:4000/build/nami.js");
  const simulationParameters = {
    skip: 1
  }
  
  const start = (canvas) => {
    const ds = 15/60;
    const extent = {
      xmin: -121,
      xmax: 25,
      ymin: -60,
      ymax: 85,
    };
    
    const scenario = {
      ...extent,
      waveWidth: parseInt((extent.xmax-extent.xmin)/ds),
      waveHeight: parseInt((extent.ymax-extent.ymin)/ds),
      earthquake: {
        ce: -12, //centroid N coordinate, e
        cn: -35,
        Mw: 9.0
      },
    };


    const nami = new Nami(scenario, { canvas }, {
      modelStepDidFinish: (model, thisController) => {
        if (model.discretization.stepNumber % simulationParameters.skip === 0) {
          return false;
        }
        return true;
      }
    });
  }

  onmessage = (e) => {
    switch(e.data.msg){
      case "start":
        start(e.data.elem);
        break;
      case "skip-toggle":
        simulationParameters.skip = simulationParameters.skip == 1 ? 50 : 1;
        postMessage(simulationParameters);
      default:
    }



  }
}

const createWorker = fn => {
  const URL = window.URL || window.webkitURL;
  return new Worker(URL.createObjectURL(new Blob(["(" + fn + ")()"])));
};

const worker = createWorker(workerProgram);

let timeEl;
let t1 = Date.now();
let t2;
function tickTimer() {
  t2 = Date.now();
  diff = t2 - t1;
  t1 = t2;
  timeEl.innerText = `approx ${(1000/diff).toFixed(1)} fps`;
  window.timerRafId = window.requestAnimationFrame(tickTimer);
}

const mainProgram = () => {
  const canvas = document.querySelector("#simulation-canvas");
  canvas.width = canvas.offsetWidth;
  canvas.height = canvas.offsetHeight;

  const offscreen = canvas.transferControlToOffscreen();
  worker.postMessage({ msg: "start", elem: offscreen }, [offscreen]);

  timeEl = document.getElementById('timer');

  requestAnimationFrame(tickTimer);
  
}

const cycleSpeed = () => {
  worker.postMessage({msg: "skip-toggle"})
}

worker.addEventListener("message", (message) => {
  document.querySelector("#speed-counter").innerText = `Skipping: ${message.data.skip-1} frames`;
})