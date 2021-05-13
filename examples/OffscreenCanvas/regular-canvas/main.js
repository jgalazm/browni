
const workerProgram = (canvas) => {
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

const mainProgram = () => {
  const canvas = document.querySelector("#simulation-canvas");
  canvas.width = canvas.offsetWidth;
  canvas.height = canvas.offsetHeight;

  workerProgram(canvas);
}


