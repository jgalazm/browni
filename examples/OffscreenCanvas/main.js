
const workerProgram = () => {
  importScripts("http://localhost:4000/build/nami.js");
  onmessage = (e) => {

    const canvas = e.data.elem;
    console.log(canvas);
    const scenario = {
      xmin: -121,
      xmax: 0,
      ymin: -60,
      ymax: 60,
      earthquake: {
        ce: -12, //centroid N coordinate, e
        cn: -35,
        Mw: 9.0
      },
    };


    const nami = new Nami(scenario, {canvas}, {
      modelSimulationWillStart: (model, thisController) => {
        thisController.paused = true;
    }
    });

  }
}


const createWorker = fn => {
  const URL = window.URL || window.webkitURL;
  return new Worker(URL.createObjectURL(new Blob(["(" + fn + ")()"])));
};

const mainProgram = () => {
  const canvas = document.querySelector("#simulation-canvas");
  canvas.width = canvas.offsetWidth;
  canvas.height = canvas.offsetHeight;

  worker = createWorker(workerProgram);
  const offscreen = canvas.transferControlToOffscreen();
  worker.postMessage({ msg: "start", elem: offscreen }, [offscreen]);
}


