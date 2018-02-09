TsunamiLab.prototype.Model = function (modelParameters) {
    
    // user parameters
    let bathymetry = modelParameters.bathymetry,
        finiteFaultParameters = modelParameters.finiteFaultParameters,
        numericalMeshResolution = modelParameters.numericalMeshResolution,
        numericalTimeStep = modelParameters.numericalTimeStep,
        CFL = modelParameters.CFL,
        coriolisFlag = modelParameters.coriolisFlag,
        pcolorMapResolution = modelParameters.pcolorMapResolution,
        colormap = modelParameters.colormap,
        pois = modelParameters.pois;
    
    // constants
    
    let R_earth = 6378000.0,
        rad_deg = 0.01745329252,
        rad_min = 0.000290888208665721,
        cori_w = 7.2722e-5;
    var shaders = shadersCode('tsunami');

    // other (global) variables
    
    let poisTime = [];
    let gridValues = [];
    let gridsTime  = []
    let simulationData = {};
    let time = 0.0;

    //---------------------------------------------
    // simulation scene
    //---------------------------------------------

    // renderer

    var rendererArguments = {
        alpha: true,
        preserveDrawingBuffer: true
    }

    var renderer = new THREE.WebGLRenderer(rendererArguments);
    renderer.setClearColor(0xffffff, 0.0);
    renderer.setSize(numericalMeshResolution.width, numericalMeshResolution.height);
    
    $(renderer.domElement).css({
        'width': pcolorMapResolution.width,
        'height': pcolorMapResolution.height
    });


    var scene = new THREE.Scene();

    var camera = new THREE.OrthographicCamera(
        -0.5 * numericalMeshResolution.width, 
        0.5 * numericalMeshResolution.width,
        0.5 * numericalMeshResolution.height, 
        -0.5 * numericalMeshResolution.height,
        -500, 1000);

    var simulation = { // TODO: check this object ......
        toggleBuffer1: false,
        toggleMaxHeightsBuffer1: false,
        paused: true,

        uniforms: {
            texel: {
                type: "v2",
                value: undefined
            },
            delta: {
                type: "v2",
                value: undefined
            },
            colormap: {
                type: "v4v",
                value: colormap
            },
            texel: { type: "v2", value: undefined },
            delta: { type: "v2", value: undefined },
            tSource: { type: "t", value: undefined },
            tMaxHeights: { type: "t", value: undefined },
            tVisualization: { type: "t", value: undefined },
            tBati: { type: "t", value: undefined },

            //sim params
            //TODO Integrar a fragment shader las constantes
            rad_min: { type: 'f', value: rad_min },
            rad_deg: { type: 'f', value: rad_deg },
            gx: { type: 'f', value: 0.01 },
            dx: { type: 'f', value: undefined },
            dy: { type: 'f', value: undefined },
            RR: { type: 'f', value: undefined },
            RS: { type: 'f', value: undefined },
            g: { type: 'f', value: 9.81 },
            xmin: { type: "f", value: 0.0 },
            xmax: { type: "f", value: 1.0 },
            ymin: { type: "f", value: 0.0 },
            ymax: { type: "f", value: 1.0 },
            zmin: { type: "f", value: 0.0 },
            zmax: { type: "f", value: 1.0 },

            //fault params
            L: { type: 'f', value: undefined },
            W: { type: 'f', value: undefined },
            depth: { type: 'f', value: undefined },
            slip: { type: 'f', value: undefined },
            strike: { type: 'f', value: undefined },
            dip: { type: 'f', value: undefined },
            rake: { type: 'f', value: undefined },
            U3: { type: 'f', value: undefined },
            cn: { type: 'f', value: undefined },   //centroid N coordinate, 18zone
            ce: { type: 'f', value: undefined },    //centroid E coordinate, 18zone

            //misc
            pause: { type: 'i', value: 0 }
        }

    }

    var simulationTextureBuffer1 = new THREE.WebGLRenderTarget(
        numericalMeshResolution.width, numericalMeshResolution.height,
        {
            minFilter: THREE.LinearFilter,
            magFilter: THREE.LinearFilter,
            format: THREE.RGBAFormat,
            type: THREE.FloatType,
            wrapS: THREE.RepeatWrapping,
            wrapT: THREE.ClampToEdgeWrapping,
            needsUpdate: true
        });

    var simulationTextureBuffer2 = new THREE.WebGLRenderTarget(
        numericalMeshResolution.width, numericalMeshResolution.height,
        {
            minFilter: THREE.LinearFilter,
            magFilter: THREE.LinearFilter,
            format: THREE.RGBAFormat,
            type: THREE.FloatType,
            wrapS: THREE.RepeatWrapping,
            wrapT: THREE.ClampToEdgeWrapping,
            needsUpdate: true
        });

    var maxHeightsTextureBuffer1 = new THREE.WebGLRenderTarget(
        numericalMeshResolution.width, numericalMeshResolution.height,
        {
            minFilter: THREE.LinearFilter,
            magFilter: THREE.LinearFilter,
            format: THREE.RGBAFormat,
            type: THREE.FloatType,
            wrapS: THREE.RepeatWrapping,
            wrapT: THREE.ClampToEdgeWrapping,
            needsUpdate: true
        });

    var maxHeightsTextureBuffer2 = new THREE.WebGLRenderTarget(
        numericalMeshResolution.width, numericalMeshResolution.height,
        {
            minFilter: THREE.LinearFilter,
            magFilter: THREE.LinearFilter,
            format: THREE.RGBAFormat,
            type: THREE.FloatType,
            wrapS: THREE.RepeatWrapping,
            wrapT: THREE.ClampToEdgeWrapping,
            needsUpdate: true
        });

    simulation.simulationTextureBuffer1 = simulationTextureBuffer1;
    simulation.simulationTextureBuffer2 = simulationTextureBuffer2;
    simulation.maxHeightsTextureBuffer1 = maxHeightsTextureBuffer1;
    simulation.maxHeightsTextureBuffer2 = maxHeightsTextureBuffer2;
    simulation.uniforms.delta.value = new THREE.Vector2(1 / numericalMeshResolution.width, 1 / numericalMeshResolution.height);



    var materials = {
        initialMaterial: new THREE.ShaderMaterial({
            uniforms: simulation.uniforms,
            vertexShader: shaders.vshader,
            fragmentShader: shaders.iFshader,
            transparent: false
        }),

        modelMaterial: new THREE.ShaderMaterial({
            uniforms: simulation.uniforms,
            vertexShader: shaders.vshader,
            fragmentShader: shaders.mFshader,
            transparent: false
        }),

        screenMaterial: new THREE.ShaderMaterial({
            uniforms: simulation.uniforms,
            vertexShader: shaders.vshader,
            fragmentShader: shaders.sFshader,
            transparent: true
        }),
        maxHeightsMaterial: new THREE.ShaderMaterial({
            uniforms: simulation.uniforms,
            vertexShader: shaders.vshader,
            fragmentShader: shaders.maxHeightsFragmentShader
        }),
        phongMaterial: new THREE.MeshBasicMaterial({
            map: bathymetry.texture
        })
    };

    var objects = {
        planeScreen: new THREE.Mesh(
            new THREE.PlaneGeometry(numericalMeshResolution.width, numericalMeshResolution.height),
            materials.phongMaterial
        )
    };

    scene.add(camera);

    scene.add(objects.planeScreen);


    //---------------------------------------------
    // renderings
    //---------------------------------------------


    var renderSimulation = function () {

        objects.planeScreen.material = materials.modelMaterial;
        
            time = time + simulationData.timeStep;

            if (simulation.toggleBuffer1) {
                simulation.uniforms.tSource.value =
                    simulation.simulationTextureBuffer1.texture;

                renderer.render(
                    scene,
                    camera,
                    simulation.simulationTextureBuffer2,
                    true
                );
            }
            else {
                simulation.uniforms.tSource.value =
                    simulation.simulationTextureBuffer2.texture;

                renderer.render(
                    scene,
                    camera,
                    simulation.simulationTextureBuffer1,
                    true
                );
            }

            simulation.toggleBuffer1 = !simulation.toggleBuffer1;
    };

    var renderMaxHeights = function () {

        objects.planeScreen.material = materials.maxHeightsMaterial;

        if (simulation.toggleMaxHeightsBuffer1) {
            simulation.uniforms.tMaxHeights.value =
                simulation.maxHeightsTextureBuffer1.texture;

            renderer.render(
                scene,
                camera,
                simulation.maxHeightsTextureBuffer2,
                true
            );
        }
        else {
            simulation.uniforms.tMaxHeights.value =
                simulation.maxHeightsTextureBuffer2.texture;

            renderer.render(
                scene,
                camera,
                simulation.maxHeightsTextureBuffer1,
                true
            );
        }

        simulation.toggleMaxHeightsBuffer1 = !simulation.toggleMaxHeightsBuffer1;
    };

    var renderScreenVoid = function () {
        // TODO: change "screen" for pcolorMap
        objects.planeScreen.material = materials.screenMaterial;
        simulation.uniforms.tVisualization.value =
            simulation.uniforms.tSource.value;

        renderer.render(
            scene,
            camera
        );
    }

    var renderMaxHeightsVisualization = function () {
        // TODO: change visualization to pcolorMap
        objects.planeScreen.material = materials.screenMaterial;
        simulation.uniforms.tVisualization.value =
            simulation.uniforms.tMaxHeights.value;

        renderer.render(
            scene,
            camera
        );
    }

    var renderScreen = function () {
        // TODO: change Screen to pcolorMap
        // TODO: change render to exportPNG
        renderScreenVoid();
        return renderer.domElement.toDataURL("image/png", 0.98);
    };


    //---------------------------------------------
    // setters and getters
    //---------------------------------------------

    var setFiniteFaultParameters = function (parameters) {
        finiteFaultParameters = parameters;
    }

    var setInitialCondition = function (faultParameters) {
        time = 0;

        setFiniteFaultParameters(faultParameters);

        simulation.uniforms.L.value = faultParameters.L;
        simulation.uniforms.W.value = faultParameters.W;
        simulation.uniforms.depth.value = faultParameters.depth;
        simulation.uniforms.slip.value = faultParameters.slip;
        simulation.uniforms.strike.value = faultParameters.strike;
        simulation.uniforms.dip.value = faultParameters.dip;
        simulation.uniforms.rake.value = faultParameters.rake;
        simulation.uniforms.U3.value = faultParameters.U3;
        simulation.uniforms.cn.value = faultParameters.cn;
        simulation.uniforms.ce.value = faultParameters.ce;

        objects.planeScreen.material = materials.initialMaterial;
        renderer.render(
            scene,
            camera,
            simulation.simulationTextureBuffer1,
            true
        );
        renderer.render(
            scene,
            camera,
            simulation.simulationTextureBuffer2,
            true
        );
        simulation.uniforms.tSource.value = simulation.simulationTextureBuffer1.texture;
        renderScreenVoid();

    }

    let setSimulation = function (finitefaultParameters) {
        // bathymetry

        simulation.uniforms.xmin.value = bathymetry.metadata.xyzmin[0];
        simulation.uniforms.ymin.value = bathymetry.metadata.xyzmin[1];
        simulation.uniforms.zmin.value = bathymetry.metadata.xyzmin[2];
        simulation.uniforms.xmax.value = bathymetry.metadata.xyzmax[0];
        simulation.uniforms.ymax.value = bathymetry.metadata.xyzmax[1];
        simulation.uniforms.zmax.value = bathymetry.metadata.xyzmax[2];

        if (simulation.uniforms.zmin.value > 0.0) {
            var e = new Error('zmin positive everywhere on bathymetry file');
            throw e;
        }

        // numerical grid / domain

        var simNx = numericalMeshResolution.width;
        var simNy = numericalMeshResolution.height;

        planeHeight = 1.0;
        planeWidth = planeHeight * simNx / simNy;
        simulation.uniforms.texel.value = new THREE.Vector2(1 / simNx, 1 / simNy)

        var dx = (simulation.uniforms.xmax.value - simulation.uniforms.xmin.value) / (simNx) * 60.0;
        var dy = (simulation.uniforms.ymax.value - simulation.uniforms.ymin.value) / (simNy) * 60.0;
        simulation.uniforms.dx.value = dx;
        simulation.uniforms.dy.value = dy;

        var lat_max = 85;//Math.max(Math.abs(ymin),Math.abs(ymax));
        var dx_real = R_earth * Math.cos(lat_max * rad_deg) * dx * rad_min;
        var dy_real = R_earth * dy * rad_min;

        var dt = 0.5 * Math.min(dx_real, dy_real) / Math.sqrt(-9.81 * simulation.uniforms.zmin.value);

        // shallow water equation constants
        simulation.uniforms.RR.value = dt / R_earth;
        simulation.uniforms.RS.value = 9.81 * simulation.uniforms.RR.value;

        simulation.uniforms.tBati.value = bathymetry.texture;


        simulationData = { 
            // TODO: should this be so deep in the model?
            // tsunamiLab.model.simulationData.gridSize 
            // or tsunamiLab.model.gridSize ??
            get gridSize(){
                return [simNx, simNy]
            },
            get bbox(){
                return [
                    [simulation.uniforms.xmin.value,simulation.uniforms.ymin.value],
                    [simulation.uniforms.xmax.value,simulation.uniforms.ymax.value]
                ]
            },
            get CFL(){
                var dt = simulation.uniforms.RR.value * R_earth;
                var celerity = Math.sqrt(-9.81 * simulation.uniforms.zmin.value);
                var cellSize =  Math.min(dx_real, dy_real);
                return dt*celerity/cellSize;
            },
            set CFL(newCFL){
                
               if(newCFL>1){

                    console.warn('CFL > 1 may induce numerical instabilities');

                }

                var celerity = Math.sqrt(-9.81 * simulation.uniforms.zmin.value);
                
                var cellSize =  Math.min(dx_real, dy_real);

                var newTimeStep = newCFL*cellSize/celerity;

                this.timeStep = newTimeStep;

            },
            get timeStep(){
                /*
                    Returns current used timestep (in seconds).
                */
                return simulation.uniforms.RR.value * R_earth;
            },
            set timeStep(timeStep){
                /*
                    Sets "timeStep" (seconds) as the current timeStep
                */

                

                simulation.uniforms.RR.value = timeStep/R_earth;

                simulation.uniforms.RS.value = 9.81 * simulation.uniforms.RR.value;

                var cfl = this.CFL;

                if(cfl>1.0){
                    
                    console.warn('CFL > 1 with new timeStep = ',timeStep,
                    'may induce numerical instabilities');

                }
            },

            get currentIterationTime(){
                /* 
                Returns the time (in minutes) of the simulation at the current iteration
                */

                return time/60.0;
            
            },

            get nextIterationTime(){
                /*
                Returns the time (in minutes)  of the simulation at the next iteration
                */

                return (time+this.timeStep )/ 60.0;
            }

        };

        setInitialCondition(finiteFaultParameters);
    }

    var getSimulationPixels = function (left, top, width, height) {
        var gridSize = simulationData.gridSize;
        var pixelData = new Float32Array(width * height * 4);
        var framebuffer = renderer.properties.get(simulation.simulationTextureBuffer1);
        framebuffer = framebuffer.__webglFramebuffer;

        var gl = getCanvas().getContext('webgl');
        gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
        gl.readPixels(left, top, width, height, gl.RGBA, gl.FLOAT, pixelData);
        return pixelData;
    }
    
    var createTextureFromData = function ( width, height, data ) {

        /* Receives an array of size width*height*4 in GL texture format
        and returns a THREE.js float texture object with this data as source.
        
        Example usage:
        
        var gl = tsunamiLab.model.getCanvas().getContext('webgl');
        data = [1.1, 1.2, 1.3, 1.4, 
            2.1, 2.2, 2.3, 2.4, 
            3.1, 3.2, 3.3, 3.4, 
            4.1, 4.2, 4.3, 4];
        t = createTextureFromData ( 2, 2, data )

        // convert to frame buffer
        fb = gl.createFramebuffer();  
        gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
        gl.framebufferTexture2D(
            gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0,
            gl.TEXTURE_2D, t.__webglTexture, 0);
    
        //read the framebuffer and unbind
        pixelData = new Float32Array(2 * 2 * 4);
        gl.readPixels(0, 0, 2, 2, gl.RGBA, gl.FLOAT, pixelData);
        gl.bindFramebuffer(gl.FRAMEBUFFER, null); 
        
        // expected output:
        // [3.0999999046325684, 3.200000047683716, 3.299999952316284, 3.4000000953674316, 
        //  4.099999904632568, 4.199999809265137, 4.300000190734863, 4, 
        //  1.100000023841858, 1.2000000476837158, 1.2999999523162842, 1.399999976158142, 
        //  2.0999999046325684, 2.200000047683716, 2.299999952316284, 2.4000000953674316]
        
        // get the numerical error
        pixelData.forEach(function(val,ind){
            c = ind%4;
            n = Math.floor(ind/4);
            console.log(pixelData[4*n+c], data[4*((n+2)%4)+c],pixelData[4*n+c]- data[4*((n+2)%4)+c]);
        })

        // about 10^-7

        */ 

        if (!gl.getExtension("OES_texture_float")) {
           throw("Requires OES_texture_float extension");
        }
        texture = new THREE.Texture( );
        texture.needsUpdate = false;
        texture.__webglTexture = gl.createTexture();

        gl.bindTexture( gl.TEXTURE_2D, texture.__webglTexture );

        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.FLOAT, new Float32Array(data) );
        texture.__webglInit = false;

        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE );
        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE );
        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR );
        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR );
        gl.generateMipmap( gl.TEXTURE_2D );
        gl.bindTexture( gl.TEXTURE_2D, null )

        return texture;
        
     }    



    var getCanvas = function () {
        return renderer.domElement;
    }

    var getFiniteFaultParameters = function () {
        return finiteFaultParameters;
    }

    var setColormap = function (colormap) {
        simulation.uniforms.colormap.value = colormap;
    };

    
    let setPOIs = function(){
        /*
            Initializes the list of pois in the model
        */
        Object.keys(pois).forEach(function(poi){
            var dlon = simulation.uniforms.dx.value;
            var dlat = simulation.uniforms.dy.value;
            var bbox = simulationData.bbox;
            var lowerLeftCorner = bbox[0];
            var upperRightCorner = bbox[1];

            var i = Math.floor((pois[poi].location[0]-lowerLeftCorner[0])/(dlon/60.0)+0.5);
            var j = Math.floor((pois[poi].location[1]-lowerLeftCorner[1])/(dlat/60.0)+0.5);

            pois[poi].pixel = [i,j];
            pois[poi].surface = [];
        });
    }

    let storePOISValues = function(){
        
        // save current time in minutes

        poisTime.push(simulationData.currentIterationTime);

        // store water surface elevation at points of interest (POIs)

        Object.keys(pois).forEach(function(poi){
            var pixelsurface = getSimulationPixels(pois[poi].pixel[0], pois[poi].pixel[1], 1, 1)[0];

            pois[poi].surface.push(pixelsurface);
            
        });
    }

    var store2DHeightValues = function(){

        // save current time in minutes

        gridsTime.push(simulationData.currentIterationTime);

        // get grid values 

        var pixels = getSimulationPixels(0, 0, numericalMeshResolution.width, numericalMeshResolution.height);

        // store them in a gridded array

        var values = [];

        for(var i=0; i< numericalMeshResolution.width; i++){

            values.push([]);

            for(var j=0; j< numericalMeshResolution.height; j++){

                var surface = pixels[4*(j*numericalMeshResolution.width + i)];

                values[i].push(surface);

            }

        }

        gridValues.push(values);

    }

    // 
    var convertGridToString = function(gridHeights){
        /*
            Takes a gridded array and returns its string encoded
            representation as csv
            
            Params.
            gridHeights: n x m array

            Returns:
            outputSting: encoded string for storing in a file

            Example:

            outputString = convertGridToString(gridHeights);
            var link = document.createElement("a");
            link.href =  outputString;
            link.download = "tlab2D00";
            link.click();

        */


        var outputString = '';
        for(var i =0; i<gridHeights.length; i++){
            for(var j=0; j<gridHeights[0].length; j++){
                outputString += gridHeights[i][j] + ' ';
            }
            outputString += '\n';
        }
        
        var outputData = new Blob([outputString], { type: 'text/csv' }); 
        var outputURL = URL.createObjectURL(outputData);
        

        return outputURL;
    }


    setColormap(colormap);
    setSimulation(finiteFaultParameters);
    setPOIs();
    storePOISValues();
    store2DHeightValues();

    return {
        renderSimulation,
        setSimulation,
        setInitialCondition,
        renderScreen,
        renderScreenVoid,
        simulation,
        zlimits: [simulation.uniforms.zmin.value, simulation.uniforms.zmax.value],
        simulationData,
        renderer,
        getCanvas,
        getFiniteFaultParameters,
        getSimulationPixels,
        renderMaxHeights,
        renderMaxHeightsVisualization,
        storePOISValues,
        pois,
        poisTime,
        store2DHeightValues,
        gridValues,
        numericalMeshResolution,
        convertGridToString,
        gridsTime
    };
}