
const workerProgram = () => {
  importScripts("http://localhost:4000/build/nami.js");
  onmessage = (e) => {

    const canvas = e.data.elem;
    console.log(canvas);
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


    const nami = new Nami(scenario, {canvas});

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


