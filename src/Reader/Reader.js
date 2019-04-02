import { getArrayFromImage, getArrayFromFile, getStringFromFile } from "./FileUtils";

const Reader = function(data, outputData) {
  let loadBathymetry = function(resolve, bathymetryInput, bathymetryMetadata, binaryBathymetry) {
    if (typeof bathymetryInput == "object") {
      // assume bathymetry is an array
      resolve(bathymetryInput);
    } else if (
      bathymetryInput.slice(-3) === "png" ||
      bathymetryInput.slice(-3) === "jpg"
    ) {
      if (!bathymetryMetadata) {
        throw new Error(
          "Must define data.bathymetryMetadata when using image format bathymetry"
        );
      }
      fetch(bathymetryInput)
      .then(r => r.blob())
      .then(blob => {
        const img = URL.createObjectURL(blob);
        const imgElement = document.createElement('img');
        imgElement.setAttribute('src', img);
        imgElement.onload = () =>{
          resolve(getArrayFromImage(imgElement, bathymetryMetadata));
        };
      });

    } else {
      getArrayFromFile(
        bathymetryInput,
        array => {
          resolve(array);
        },
        binaryBathymetry ? "binary" : "ascii"
      );
    }
  };

  let loadInitialCondition = function(resolve) {
    /* Detects if initialSurface, earthquake or asteroid is provided, assuming
        the user knows the right format.
        Otherwise throws an error */

    if (data.initialSurface != undefined) {
      getArrayFromFile(
        data.initialSurface.file,
        function(array) {
          data.initialSurface.array = array;

          initialSurfaceReady = true;

          resolve(data.initialSurface.array);
        },
        "ascii"
      );
    } else if (data.earthquake != undefined) {
      if (typeof data.earthquake === "string") {
        getStringFromFile(data.earthquake, fileString => {
          let earthquake = fileString
            .replace(String.fromCharCode(13), "")
            .split("\n");
          let keys = earthquake[0].split(",");
          let key2column = {};
          for (let i = 0; i < keys.length; i++) {
            key2column[keys[i]] = i;
          }

          earthquake.shift();

          earthquake = earthquake.filter(val => val.length > 0);

          for (let i = 0; i < earthquake.length; i++) {
            let finiteFault = earthquake[i].split(",");
            let earthquakeDict = {};
            keys.map(key => {
              earthquakeDict[key] =
                key != "reference"
                  ? parseFloat(finiteFault[key2column[key]])
                  : finiteFault[key2column[key]];
            });
            earthquake[i] = earthquakeDict;
          }

          data.earthquake = earthquake;

          resolve(earthquake);
        });
      } else {
        resolve(data.earthquake);
      }
    } else if (data.asteroid !== undefined) {
      resolve(data.asteroid);
    } else {
      throw "Need a valid data.initialSurface or data.finiteFault as input";
    }
  };

  // domain
  let domain = {
    equations: data.equations,
    coordinates: data.coordinates,
    xmin: data.xmin,
    xmax: data.xmax,
    ymin: data.ymin,
    ymax: data.ymax,
    isPeriodic: data.isPeriodic !== undefined ? data.isPeriodic : 0
  };

  // discretization
  let [nx,ny] = [-1, -1];
  if(data.gridSize !== undefined){
    nx = (domain.xmax-domain.xmin)/gridSize + 1;
    ny = (domain.ymax-domain.ymin)/gridSize + 1;
  }
  else if( data.waveWidth !== undefined && data.waveHeight !== undefined){
    [nx,ny] = [data.waveWidth, data.waveHeight];
  }
  else if(data.bathymetry === undefined){
    // 20 minutes resolution, 1080 x 480 grid for the whole globe
    nx = (domain.xmax-domain.xmin)*3 + 1;
    ny = (domain.ymax-domain.ymin)*3 + 1;
  }
  else{
    throw "Cannot determine grid size from input";
  }

  let discretization = {
    numberOfCells: [nx, ny],
    dt: undefined,
    stepNumber: 0
  };

  if (domain.coordinates === undefined){
    domain.coordinates = "spherical";
  }

  if (domain.coordinates == "cartesian") {
    discretization.dx =
      (domain.xmax - domain.xmin) / (discretization.numberOfCells[0] - 1);
    discretization.dy =
      (domain.ymax - domain.ymin) / (discretization.numberOfCells[1] - 1);
  } else if (domain.coordinates == "spherical") {
    domain.xmin = domain.xmin;
    domain.xmax = domain.xmax;
    discretization.dlon =
      (60 * (domain.xmax - domain.xmin)) /
      (discretization.numberOfCells[0] - 1);
    discretization.dlat =
      (60 * (domain.ymax - domain.ymin)) /
      (discretization.numberOfCells[1] - 1);
  }
  isPeriodic
  if (domain.equations === undefined) {
    domain.equations = "linear";
  }

  // slab

  let slab = new Promise( (resolve, reject) => {
    if(data.slab) {
      resolve(data.slab);
      return;
    }

    fetch("https://unpkg.com/namijs@1.0.3/assets/SlabModels.json").then( response => {
      resolve(response.json());
    });
  });

  // bathymetry
  let bathymetry = {
    array: undefined,
    image: undefined,
    textureId: 0,
    texture: undefined, // to be loaded at start()
    extent: {
      xmin: domain.xmin,
      xmax: domain.xmax,
      ymin: domain.ymin,
      ymax: domain.ymax
    }
  };

  if (data.bathymetryExtent) {
    bathymetry.extent = {
      xmin: data.bathymetryExtent.xmin,
      xmax: data.bathymetryExtent.xmax,
      ymin: data.bathymetryExtent.ymin,
      ymax: data.bathymetryExtent.ymax
    };
  }

  if (
    domain.xmin < bathymetry.extent.xmin ||
    domain.ymin < bathymetry.extent.ymin ||
    domain.xmax > bathymetry.extent.xmax ||
    domain.ymax > bathymetry.extent.ymax
  ) {
    throw `Bathymetry extent should include domain extent 
           domain:${domain}, bathymetry.extent:${bathymetry.extent}`;
  }
  
  let bathymetryInput = data.bathymetry; 
  let bathymetryMetadata = data.bathymetryMetadata;
  if( data.bathymetry === undefined){
    bathymetryInput = 'https://unpkg.com/namijs@1.0.3/assets/bathymetry.png';
    bathymetryMetadata = {
      zmin: -6709,
      zmax: 10684
    };
  }

  bathymetry.array = new Promise((resolve, reject) => {
    loadBathymetry(resolve, bathymetryInput, bathymetryMetadata, data.binaryBathymetry);
  });

  // initial condition

  let initialCondition = new Promise((resolve, reject) => {
    loadInitialCondition(resolve);
  });

  const pcolorDisplay = {
    width: outputData.displayWidth ? outputData.displayWidth : discretization.numberOfCells[0],
    height: outputData.displayHeight ? outputData.displayHeight : discretization.numberOfCells[1]
  };

  const displayOption = outputData.displayOption
    ? outputData.displayOption
    : "heights";

  const pois = output.pois ? output.pois : {};

  const cmax = 0.1;
  const cmin = -0.1;
  let defaultColormap = {
    thresholds: [
      0.0 * (cmax - cmin) + cmin,
      0.06666666666666667 * (cmax - cmin) + cmin,
      0.13333333333333333 * (cmax - cmin) + cmin,
      0.2 * (cmax - cmin) + cmin,
      0.26666666666666666 * (cmax - cmin) + cmin,
      0.3333333333333333 * (cmax - cmin) + cmin,
      0.4 * (cmax - cmin) + cmin,
      0.49 * (cmax - cmin) + cmin,
      0.5 * (cmax - cmin) + cmin,
      0.51 * (cmax - cmin) + cmin,
      0.6666666666666666 * (cmax - cmin) + cmin,
      0.7333333333333333 * (cmax - cmin) + cmin,
      0.8 * (cmax - cmin) + cmin,
      0.8666666666666667 * (cmax - cmin) + cmin,
      0.9333333333333333 * (cmax - cmin) + cmin,
      1.0 * (cmax - cmin) + cmin
    ],

    rgba: [
      [0.001462, 0.000466, 0.013866, 1],
      [0.046915, 0.030324, 0.150164, 0.8],
      [0.142378, 0.046242, 0.308553, 0.8],
      [0.258234, 0.038571, 0.406485, 0.8],
      [0.366529, 0.071579, 0.431994, 0.8],
      [0.472328, 0.110547, 0.428334, 0.9],
      [0.578304, 0.148039, 0.404411, 0.8],
      [0.682656, 0.189501, 0.360757, 0.4],
      [0.780517, 0.243327, 0.299523, 0],
      [0.865006, 0.316822, 0.226055, 0.4],
      [0.929644, 0.411479, 0.145367, 0.8],
      [0.970919, 0.522853, 0.058367, 0.9],
      [0.987622, 0.64532, 0.039886, 0.8],
      [0.978806, 0.774545, 0.176037, 0.8],
      [0.950018, 0.903409, 0.380271, 0.8],
      [0.988362, 0.998364, 0.644924, 1]
    ]
  };

  let colormap =
    output.colormap !== undefined ? output.colormap : defaultColormap;

  // flatten the array
  colormap.rgba = colormap.rgba.reduce((a, b) => {
    return a.concat(b);
  });

  return {
    domain,
    discretization,
    bathymetry,
    initialCondition,
    slab,
    pcolorDisplay,
    displayOption,
    pois,
    colormap
  };
};

export default Reader;
