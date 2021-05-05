import { getLengthWidthSlip } from "./Earthquake";
import Earthquake from "./renderers/Earthquake/Earthquake";
import SphericalShallowWater from "./renderers/SphericalShallowWater/SphericalShallowWater";
import Display from "./renderers/Display/Display";
import MaxHeights from "./renderers/MaxHeights/MaxHeights";

let Model = function(data) {
  let gl, isWebGL2;

  let wave, maxHeights;
  let {
    domain,
    bathymetry,
    discretization,
    initialSurface,
    earthquake,
    asteroid,
    slab,
    pcolorDisplay,
    displayOption,
    pois,
    colormap
  } = data;

  /* 
        Start WebGL context

    */

  let canvas;
  if (data.canvas != undefined) {
    canvas = data.canvas;
  } else {
    canvas = document.createElement("canvas");
  }
  canvas.width = pcolorDisplay.width;
  canvas.height = pcolorDisplay.height;

  try {
    gl = canvas.getContext("webgl2", { premultipliedAlpha: false });
    isWebGL2 = !!gl;
    if (!isWebGL2) {
      gl =
        canvas.getContext("webgl", { premultipliedAlpha: false }) ||
        canvas.getContext("experimental-webgl", { premultipliedAlpha: false });
    }
    gl.viewportWidth = canvas.width;
    gl.viewportHeight = canvas.height;
  } catch (e) {}
  if (!gl) {
    alert("Could not initialise Webgl.");
  }

  gl.enable(gl.DEPTH_TEST);

  gl.clearColor(0.0, 0.0, 0.0, 1.0); //default color for any fbo, canvas included

  const modelState = {
    discretization,
    bathymetry,
    domain
  };

  const earthquakeModel = new Earthquake(gl);
  const sphericalShallowWaterModel = new SphericalShallowWater(gl);
  const displayStep = new Display(gl);
  const maxHeightsStep = new MaxHeights(gl);
  /*
        WebGL tools

    */

  let createBuffers = function() {
    let vertexPositionBufferHandle = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexPositionBufferHandle);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      new Float32Array([-1, -1, -1, 1, 1, 1, 1, -1]),
      gl.STATIC_DRAW
    );

    let facesBufferHandle = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, facesBufferHandle);
    gl.bufferData(
      gl.ELEMENT_ARRAY_BUFFER,
      new Uint16Array([0, 1, 2, 0, 2, 3]),
      gl.STATIC_DRAW
    );

    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
  };

  let createTextureFromData = function(
    width,
    height,
    data,
    textureId,
    internalFormat,
    format,
    type
  ) {
    // creates a texture from an array "data" of size "width"*height"*4 with the given webgl formats and id

    if (
      !gl.getExtension("OES_texture_float") &&
      !gl.getExtension("EXT_color_buffer_float")
    ) {
      throw "Requires OES_texture_float extension";
    }

    var texture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0 + textureId);
    gl.bindTexture(gl.TEXTURE_2D, texture);

    // 32F: diferencia con webgl/webgl2

    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      internalFormat,
      width,
      height,
      0,
      format,
      type,
      new Float32Array(data)
    );

    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    return { texture, textureId };
  };

  let createTextureFromMatrix = function(matrix, textureId) {
    let internalFormat = isWebGL2 ? gl.RGBA32F : gl.RGBA;
    let format = gl.RGBA;
    let type = gl.FLOAT;

    let raveledMatrix = [];
    for (let j = 0; j < matrix.length; j++) {
      for (let i = 0; i < matrix[0].length; i++) {
        raveledMatrix.push(matrix[j][i]);
        raveledMatrix.push(0);
        raveledMatrix.push(0);
        raveledMatrix.push(1);
      }
    }

    let texture = createTextureFromData(
      matrix[1].length,
      matrix.length,
      raveledMatrix,
      textureId,
      internalFormat,
      format,
      type
    );

    return texture;
  };

  let createTextureFromLongArray = function(longArray, textureId) {
    let internalFormat = isWebGL2 ? gl.R32F : gl.RGBA;
    let format = gl.RED;
    let type = gl.FLOAT;

    let texture = createTextureFromData(
      longArray[1],
      longArray[0],
      longArray.slice(2, longArray.length ),
      textureId,
      internalFormat,
      format,
      type
    );

    return texture;
  };

  let createFBO = function(
    textureId,
    w,
    h,
    internalFormat,
    format,
    type,
    param
  ) {
    /* 
            textureId: integer that identifies this texture such that texId+gl.TEXTURE0 is the texture unit bound to it
            w: width (pixels)
            h: height (pixels)
            internalFormat: (webgl man) Specifies the internal format of the texture.  Must be one of the following symbolic constants: GL_ALPHA, GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA.
            format: Specifies the format of the texel data. Must match internalformat. The following symbolic values are accepted: GL_ALPHA, GL_RGB, GL_RGBA, GL_LUMINANCE, and GL_LUMINANCE_ALPHA.
            type: Specifies the data type of the texel data. The following symbolic values are accepted: GL_UNSIGNED_BYTE, GL_UNSIGNED_SHORT_5_6_5, GL_UNSIGNED_SHORT_4_4_4_4, and GL_UNSIGNED_SHORT_5_5_5_1.
            (see: webgl 1.0: https://www.khronos.org/registry/OpenGL-Refpages/es2.0/xhtml/glTexImage2D.xml
                webgl 2.0: https://www.khronos.org/registry/OpenGL-Refpages/es3.0/html/glTexImage2D.xhtml)
            param: type of min/mag filter: LINEAR, NEAREST, etc.

        */
    // var internalFormat = gl.RGBA32F, format = gl.RGBA, type = gl.FLOAT, param = gl.NEAREST;

    gl.activeTexture(gl.TEXTURE0 + textureId);
    let texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, param);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, param);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      internalFormat,
      w,
      h,
      0,
      format,
      type,
      null
    ); //data will be rendered later

    let fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(
      gl.FRAMEBUFFER,
      gl.COLOR_ATTACHMENT0,
      gl.TEXTURE_2D,
      texture,
      0
    );
    gl.viewport(0, 0, w, h);
    gl.clear(gl.COLOR_BUFFER_BIT);

    let clearBuffer = () => {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
      gl.clearColor(0, 0, 0, 1);
      gl.clear(gl.COLOR_BUFFER_BIT);
    };
    return { texture, fbo, textureId, clearBuffer };
  };

  let createDoubleFBO = function(
    textureId1,
    textureId2,
    w,
    h,
    internalFormat,
    format,
    type,
    param
  ) {
    let fbo1 = createFBO(textureId1, w, h, internalFormat, format, type, param);
    let fbo2 = createFBO(textureId2, w, h, internalFormat, format, type, param);

    return {
      get first() {
        return fbo1;
      },
      get second() {
        return fbo2;
      },
      swap: function() {
        let temp = fbo1;
        fbo1 = fbo2;
        fbo2 = temp;
      },
      clearBuffers: () => {
        fbo1.clearBuffer();
        fbo2.clearBuffer();
      }
    };
  };

  let readFBOPixels = function(frameBufferObject, left, top, width, height) {
    var pixelData = new Float32Array(width * height * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBufferObject);
    gl.readPixels(left, top, width, height, gl.RGBA, gl.FLOAT, pixelData);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

    return pixelData;
  };

  let exportBuffer = function(
    fbo,
    variableIndex = 0,
    iStart = 0,
    jStart = 0,
    Lx = discretization.numberOfCells[0],
    Ly = discretization.numberOfCells[1]
  ) {
    let array = readFBOPixels(fbo, iStart, jStart, Lx, Ly);
    array = array.filter((elem, index) => {
      return (index - variableIndex) % 4 == 0;
    });

    return array;
  };

  /* Additional programs that use the results of simulations */

  let setPOIs = function() {
    /*
            Initializes the list of pois in the model
        */

    if (Object.keys(pois).length === 0) return;

    // let bathymetryTemp = exportBuffer(
    //   wave.first.fbo,
    //   3,
    //   0,
    //   0,
    //   data.waveWidth,
    //   data.waveHeight
    // );
    // bathymetryTemp = [...bathymetryTemp];

    // let bathymetry = [];

    // while (bathymetryTemp.length > 0)
    //   bathymetry.push(bathymetryTemp.splice(0, data.waveWidth));

    Object.keys(pois).forEach(function(poi) {
      const dlon = discretization.dlon;
      const dlat = discretization.dlat;
      const lowerLeftCorner = [domain.xmin, domain.ymin];

      pois[poi].location[0] = pois[poi].location[0];

      const i = Math.floor(
        (pois[poi].location[0] - lowerLeftCorner[0]) / (dlon / 60.0) + 0.5
      );
      const j = Math.floor(
        (pois[poi].location[1] - lowerLeftCorner[1]) / (dlat / 60.0) + 0.5
      );

      pois[poi].pixel = [i, j];
      pois[poi].surface = [];
      pois[poi].time = [];

      // if depth is provided and is shallow then use it, otherwise get it from the matrix
      // the texture is read in [j][i] order
      // pois[poi].depth =
      //   pois[poi].depth && pois[poi].depth < 100
      //     ? pois[poi].depth
      //     : bathymetry[j][i];

      pois[poi].shallowCorrectionFactor = 1;
      // pois[poi].closestDeepPoint = [];
      // pois[poi].closestDeepPointDepth = undefined;

      // if point is in shallow water look for closest point in deep water:
      // let closestDeepPointDistance = Infinity;
      // if (pois[poi].depth <= 100) {
      //   console.log("shallow poi:", poi);
      //   for (let j0 = 0; j0 < bathymetry.length; j0++) {
      //     for (let i0 = 0; i0 < bathymetry[0].length; i0++) {
      //       const distance = (i0 - i) * (i0 - i) + (j0 - j) * (j0 - j); // assumes uniform cartesian grid
      //       if (
      //         bathymetry[j0][i0] > 100.0 &&
      //         distance < closestDeepPointDistance
      //       ) {
      //         pois[poi].closestDeepPoint = [i0, j0];
      //         pois[poi].closestDeepPointDepth = bathymetry[j0][i0];
      //         closestDeepPointDistance = distance;
      //       }
      //     }
      //   }

      //   const d = Math.max(pois[poi].depth, 1.0);
      //   const d0 = pois[poi].closestDeepPointDepth;
      //   pois[poi].originalPixel = pois[poi].pixel;
      //   pois[poi].pixel = pois[poi].closestDeepPoint;
      //   pois[poi].shallowCorrectionFactor = Math.pow(d0 / d, 0.25);
      // }
    });
  };

  let storePOISValues = function() {
    // save current time in minutes

    // poisTime.push(simulationData.currentIterationTime);

    // store water surface elevation at points of interest (POIs)

    Object.keys(pois).forEach(function(poi) {
      var pixelSurface = readFBOPixels(
        wave.first.fbo,
        pois[poi].pixel[0],
        pois[poi].pixel[1],
        1,
        1
      )[0];
      pixelSurface = pixelSurface * pois[poi].shallowCorrectionFactor;
      pois[poi].surface.push(pixelSurface);
      pois[poi].time.push(discretization.dt * discretization.stepNumber);
    });
  };

  let runSimulationStep = function() {
    discretization.stepNumber++;

    if (domain.equations == "dispersive") {
      // choose dispersive solver
      if (domain.coordinates == "cartesian") {
        renderCartesianDispersiveProgram();
      } else if (domain.coordinates == "spherical") {
        renderDispersiveWaveEquation();
      }
    } else {
      // choose non-dispersive solver
      if (domain.coordinates == "cartesian") {
        renderCartesianProgram();
      } else if (domain.coordinates == "spherical") {
        sphericalShallowWaterModel.render(wave, modelState);
      }
    }

    maxHeightsStep.render(maxHeights, wave, modelState);

    storePOISValues();
  };

  let setTimeStep = function(options) {
    let hmax = Array.isArray(bathymetry.array[0]) ? Math.max.apply(
      null,
      bathymetry.array.map(row => {
        return Math.max.apply(null, row);
      })
    ) : bathymetry.array.reduce((prev,curr) => Math.max(prev, curr));

    let cfl;
    if (options.timeStep !== undefined) {
      discretization.dt = options.timeStep;
    } else {
      /* If dt is not given use given or predefined cfl */

      cfl = options.cfl === undefined ? 0.5 : options.cfl;
      if (domain.coordinates == "cartesian") {
        discretization.dt =
          (cfl * Math.min(discretization.dx, discretization.dy)) /
          Math.sqrt(hmax * 9.81);
      } else if (domain.coordinates == "spherical") {
        // let Rearth = 6378000.0;
        let Rearth = 6378000;
        let radPerDeg = 0.01745329252;
        let radPerMin = 0.000290888208665721;
        var latMax = Math.max(Math.abs(domain.ymax), Math.abs(domain.ymin)); //Math.max(Math.abs(ymin),Math.abs(ymax));
        var dxReal =
          Rearth *
          Math.cos(latMax * radPerDeg) *
          discretization.dlon *
          radPerMin;
        var dyReal = Rearth * discretization.dlat * radPerMin;

        var Dx = Rearth * discretization.dlon * radPerMin;
        var Dy = Rearth * discretization.dlat * radPerMin;
        var dt = (Math.cos(latMax * radPerDeg) * Dy) / Math.sqrt(9.81 * hmax);

        // let hmax = Math.max.apply(null, )
        // discretization.dt = 0.25 * Math.min(dxReal, dyReal) / Math.sqrt(9.81 * hmax);
        discretization.dt = cfl * dt;
      }
    }

    // if(discretization.dt>15)
    //     discretization.dt = 15;

    // let hmax = Math.max.apply(null, )
  };

  let getSlabParameters = (lon, lat) => {
    if (slab == undefined) {
      return;
    }

    let calcDistance = (p1, p2) => {
      return (
        (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
      );
    };
    let xslab = slab["x"][0] > 180 ? slab["x"][0] - 360 : slab["x"][0];

    let minDistance = calcDistance([xslab, slab["y"][0]], [lon, lat]);
    let minDistanceLocation = 0;
    for (let i = 1; i < slab["x"].length; i++) {
      xslab = slab["x"][i] > 180 ? slab["x"][i] - 360 : slab["x"][i];
      let newMinDistance = calcDistance([xslab, slab["y"][i]], [lon, lat]);
      if (newMinDistance < minDistance) {
        minDistance = newMinDistance;
        minDistanceLocation = i;
      }
    }
    return {
      depth: slab["depth"][minDistanceLocation],
      dip: slab["dip"][minDistanceLocation],
      strike: slab["strike"][minDistanceLocation],
      distance: minDistance,
      minDistanceLocation: minDistanceLocation
    };
  };

  let setEarthquake = () => {
    if (earthquake.length > 0) {
      for (let i = 0; i < earthquake.length; i++) {
        if (earthquake[i].lat !== undefined) {
          earthquake[i].cn = earthquake[i].lat;
        }
        if (earthquake[i].lon !== undefined) {
          earthquake[i].ce = earthquake[i].lon;
        }
        earthquake[i].U3 = 0.0;

        if (earthquake[i].rake === undefined) {
          earthquake[i].rake = 90.0;
        }

        if (earthquake[i].reference === undefined) {
          earthquake[i].reference = "center";
        }

        if (
          earthquake[i].Mw != undefined &&
          !(
            earthquake[i].L != undefined &&
            earthquake[i].W != undefined &&
            earthquake[i].slip != undefined
          )
        ) {
          const LWslip = getLengthWidthSlip(earthquake[i].Mw);
          earthquake[i].L = LWslip.L;
          earthquake[i].W = LWslip.W;
          earthquake[i].slip = LWslip.slip;

          const slabInfo = getSlabParameters(
            earthquake[i].ce,
            earthquake[i].cn
          );
          if (slabInfo) {
            earthquake[i].depth = -slabInfo.depth * 1000;
            earthquake[i].dip = slabInfo.dip;
            earthquake[i].strike = slabInfo.strike;
          }
        }
      }
    }
  };

  let initFBOs = function() {
    const internalFormat = isWebGL2 ? gl.RGBA32F : gl.RGBA;
    const format = gl.RGBA;
    const textype = gl.FLOAT;
    let textureIdHeights1 = 2;
    let textureIdHeights2 = 3;
    let textureIdMaxHeights1 = 4;
    let textureIdMaxHeights2 = 5;
    let param = gl.NEAREST;

    wave = createDoubleFBO(
      textureIdHeights1,
      textureIdHeights2,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1],
      internalFormat,
      format,
      textype,
      param
    );
    maxHeights = createDoubleFBO(
      textureIdMaxHeights1,
      textureIdMaxHeights2,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1],
      internalFormat,
      format,
      textype,
      param
    );

    if (initialSurface) {
      renderInitialProgram();
    } else if (earthquake.length > 0) {
      earthquakeModel.render(wave, modelState, earthquake);
    } else if (asteroid) {
      renderAsteroidProgram();
    }

    maxHeightsStep.render(maxHeights, wave, modelState);

    displayStep.render(wave, maxHeights, displayOption, colormap);
  };

  let start = function() {
    bathymetry.texture = createTextureFromLongArray(
      bathymetry.array,
      bathymetry.textureId
    );

    if (initialSurface != undefined) {
      initialSurface.texture = createTextureFromMatrix(
        initialSurface.array,
        initialSurface.textureId
      );

      initialSurface = Object.assign(data.initialSurface, initialSurface);
      if (
        data.initialSurface.cn === undefined ||
        data.initialSurface.ce === undefined ||
        data.initialSurface.L === undefined ||
        data.initialSurface.W === undefined
      ) {
        initialSurface.isSubrectangle = false;
      }
    }

    setTimeStep(data);

    createBuffers();

    setEarthquake();

    initFBOs();

    setPOIs();
  };

  start();

  return {
    domain,
    discretization,
    bathymetry,
    getCellBathymetry: (i, j) => {
      return exportBuffer(wave.first.fbo, 3, i, j, 1, 1);
    },
    getRectangleBathymetry: (iStart, jStart, iEnd, jEnd) => {
      let rectangle = exportBuffer(
        wave.first.fbo,
        3,
        iStart,
        jStart,
        iEnd - iStart + 1,
        jEnd - jStart + 1
      );
      rectangle = [...rectangle];
      let matrixRectangle = [];

      while (rectangle.length > 0)
        matrixRectangle.push(rectangle.splice(0, iEnd - iStart + 1));

      return matrixRectangle;
    },
    get currentTime() {
      return discretization.stepNumber * discretization.dt;
    },
    setTimeStep,
    canvas,
    get currentGridHeights() {
      return exportBuffer(wave.first.fbo);
    },
    get currentMaximumHeights() {
      return exportBuffer(maxHeights.first.fbo);
    },
    get currentArrivalTimes() {
      return exportBuffer(maxHeights.first.fbo, 1);
    },
    pois,
    runSimulationStep,
    displayPColor: () => {
      displayStep.render(wave, maxHeights, displayOption, colormap);
    },
    displayOption,
    set colors(newColors) {
      colormap.rgba = [...newColors].reduce((a, b) => {
        return a.concat(b);
      });
      let nextgrid = 0;
      gl.uniform4fv(
        displayProgram.uniforms.colormap,
        new Float32Array(colormap.rgba)
      );
    },
    get colors() {
      return colormap.rgba;
    },
    set earthquake(newEarthquake) {
      console.log(gl);
      wave.clearBuffers();
      earthquake = newEarthquake;
      setEarthquake();
      renderEarthquake();
      displayStep.render(wave, maxHeights, displayOption, colormap);
      discretization.stepNumber = 0;
    },

    get earthquake() {
      return Object.assign({}, earthquake);
    }
  };
};

export { Model };
