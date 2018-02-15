"use strict";

let NAMI = function(data, output, lifeCycle){

    let handler = this;
    this.model = {};
    this.controller = {};
    let bathymetryReady = false;
    let initialSurfaceReady = false;


    let getFile = function(url, callback){
        $.ajax({
            dataType: 'text',
            url: url,
            async: true,
            success: callback,
        });
    }
    
    let parseFile = function(data){
        data = data.split('\n');
        let arr = data.map(function(row){
            return row.split(/(\s+)/).filter( function(e) {
                return e.trim().length > 0;  
            }).map(function(val){ return parseFloat(val)});; 
        });
        arr = arr.filter(function(e){
            return e.length>0; // prevents blank lines
        })
    
        return arr;
    }

    let getArrayFromImage = function(image){
        let canvas = document.createElement('canvas');
        canvas.height = image.height;
        canvas.width = image.width;

        let ctx = canvas.getContext('2d');
        ctx.drawImage(image, 0, 0);
        let imageData = ctx.getImageData(0,0,canvas.width,canvas.height);
        
        imageData = imageData.data.filter((value, index)=>{
            return index % 4 == 0;
        });

        // convert to normal float numbers
        imageData = [...imageData];


        imageData = imageData.map((value)=>{
            return (value/255.0*(data.bathymetryMetadata.zmax-
                data.bathymetryMetadata.zmin)+
                data.bathymetryMetadata.zmin);
        });

        let matrix = [];
        while(imageData.length > 0) matrix.push(imageData.splice(0, image.width));

        return matrix;
    }
    let init = (Model) => {
        
        this.model = new this.Model(data, output);

        if (lifeCycle.dataWasLoaded !== undefined){
            lifeCycle.dataWasLoaded(this.model);
        }

        this.controller = new this.Controller(this.model,output, lifeCycle);

        this.controller.animate();
        
    }

    let loadBathymetry = function(){
        if(data.bathymetry.slice(-3)==='png' || data.bathymetry.slice(-3)==='jpg'){
            if(!data.bathymetryMetadata){
                throw new Error('Must define data.bathymetryMetadata when using image format bathymetry');
            }
            let bathymetryImage = new Image();
            bathymetryImage.onload = function(){

                data.bathymetry = {
                    array : getArrayFromImage(bathymetryImage)
                }

                bathymetryReady = true;
                
                if(initialSurfaceReady){
                    init();
                }
            }
            bathymetryImage.src = data.bathymetry;

        }
        else{
            getFile(data.bathymetry,function(fileString){
                
                data.bathymetry = {
                    array : parseFile(fileString)  
                }   
        
                bathymetryReady = true;
                
                if(initialSurfaceReady){
                    init();
                }
                
            });
        }
    }

    let loadInitialCondition = function(){
        /* Detects if initialSurface or finiteFault is provided, assuming
        the user knows the right format.
        Otherwise throws an error */
        if( data.initialSurface != undefined){

            getFile(data.initialSurface,function(fileString){
                
                data.initialSurface = {
                    array : parseFile(fileString)
                }
        
                initialSurfaceReady = true;
        
                if(bathymetryReady){           
                    init();
                }
                
            });
        }
        else if( data.earthquake != undefined){
            if(typeof(data.earthquake)==="string"){
                
                getFile(data.earthquake, fileString => {
                    let earthquake = fileString.split('\n');
                    let keys = earthquake[0].split(',');
                    let key2column = {};
                    for(let i = 0; i<keys.length; i++){
                        key2column[keys[i]] = i
                    }

                    earthquake.shift();

                    earthquake = earthquake.filter((val)=> val.length>0);

                    for(let i = 0; i<earthquake.length; i++){
                        let finiteFault = earthquake[i].split(',');
                        if(earthquake[i].length == 0){}
                        earthquake[i] ={
                            L: parseFloat(finiteFault[key2column['L']]),
                            W: parseFloat(finiteFault[key2column['W']]),
                            depth: parseFloat(finiteFault[key2column['depth']]),
                            slip: parseFloat(finiteFault[key2column['slip']]),
                            strike: parseFloat(finiteFault[key2column['strike']]),
                            dip: parseFloat(finiteFault[key2column['dip']]),
                            rake: parseFloat(finiteFault[key2column['rake']]),
                            U3: 0.0,
                            cn: parseFloat(finiteFault[key2column['lat']]),
                            ce: parseFloat(finiteFault[key2column['lon']]),
                            reference: finiteFault[key2column['reference']]
                        }
                    }
                    data.earthquake = earthquake;

                    initialSurfaceReady = true;

                    if(bathymetryReady){
                        init();
                    }

                });
            }
            else{
                
                initialSurfaceReady = true;
            }

            if(bathymetryReady && initialSurfaceReady){
                init();
            }
        }
        else {
            throw 'Need a valid data.initialSurface or data.finiteFault as input';
        }
    }

    loadBathymetry();
    loadInitialCondition();

}


NAMI.prototype.Controller = function(model,data, lifeCycle){
    let thisController = {};

    
    let controllerSimulationDidFinish, modelStepDidFinish;

    if( typeof lifeCycle !== 'undefined'){
        if( typeof lifeCycle.controllerSimulationDidFinish !== 'undefined'){
            controllerSimulationDidFinish = lifeCycle.controllerSimulationDidFinish;
        }

        if( typeof lifeCycle.modelStepDidFinish !== 'undefined'){
            modelStepDidFinish = lifeCycle.modelStepDidFinish;
        }
    }

    let downloadGridArray = function(array, time){
        
        array = array.toString();
        
        array = time.toString()+'\n'+array;

        var outputData = new Blob([array], { type: 'text/csv' }); 
        var outputURL = URL.createObjectURL(outputData);
        var link = document.createElement('a');
            
        link.href =  outputURL;
        
        link.download = 'tlab2D';

        link.click();  
    }

    let downloadPOIs = function(pois){
        var exportString = "#seconds";
        Object.keys(pois).forEach(function(poi){
            exportString +=','+poi+'(m)';
        });
        exportString+='\n';        

        
        for(let i=0;i< pois[Object.keys(pois)[0]].time.length; i++){
            exportString += pois[Object.keys(pois)[0]].time[i].toFixed(8).toString();
            
            
            Object.keys(pois).forEach(function(poi){
                exportString += ','+pois[poi].surface[i].toFixed(8).toString();
            });

            exportString += '\n';
        }
        
        exportString = 'data:text/csv;charset=utf-8,'+exportString;   
        exportString = encodeURI(exportString);     

        let link = document.createElement("a");
        link.href =  exportString;
        link.download = "tseries";
        link.click();   

        }

    let animate = () => {

        /********************* */
        model.runSimulationStep();                
            
        modelStepDidFinish(model, thisController);

        model.displayPColor();
        /********************* */


        if(model.currentTime<data.stopTime){
            
            requestAnimationFrame(animate);
            
        }
        else{

            controllerSimulationDidFinish(model, thisController);
            
        }
    
    }

    thisController = {
        animate,
        downloadCurrentGridHeights: ()=> {
            downloadGridArray(model.currentGridHeights,
                model.discretization.stepNumber*model.discretization.dt);
        },
        downloadMaximumHeights: () =>{
            downloadGridArray(model.currentMaximumHeights,
                model.discretization.stepNumber*model.discretization.dt);
        },
        downloadArrivalTimes: () =>{
            downloadGridArray(model.currentArrivalTimes,
                model.discretization.stepNumber*model.discretization.dt);
        },
        downloadAllPois: () =>{
            downloadPOIs(model.pois);
        },
        get stopTime(){
            return data.stopTime;
        }
    };

    return thisController;
}

NAMI.prototype.Model = function(data, output){
    let gl, isWebGL2;
    let vertexShader, initialShader, okadaShader, cartesianWaveShader, 
    sphericalWaveShader, maxHeightsShader,displayShader;
    
    let initialProgram, okadaProgram, cartesianWaveProgram, 
    sphericalWaveProgram, maxHeightsProgram, displayProgram  ;

    
    let domain, bathymetry, discretization, initialSurface, earthquake, wave, maxHeights, pcolorDisplay, pois, colors;
    
    let displayOption = output.displayOption ? output.displayOption : 'heights';

    // domain
    domain = {
        coordinates : data.coordinates,
        xmin : data.xmin,
        xmax : data.xmax,
        ymin : data.ymin,
      ymax : data.ymax,
      isPeriodic: data.isPeriodic !== undefined ? data.isPeriodic: 0
    };
    
    // bathymetry
    bathymetry = {
        array: data.bathymetry.array,
        image: data.bathymetry.image,
        textureId: 0,
        texture: undefined, // to be loaded at start()
    };

    // discretization
    discretization = {
        numberOfCells : [data.waveWidth, data.waveHeight],
        dt: undefined,
        stepNumber : 0
    };

    if(domain.coordinates == 'cartesian'){
        discretization.dx = (domain.xmax-domain.xmin)/(discretization.numberOfCells[0]-1)
        discretization.dy = (domain.ymax-domain.ymin)/(discretization.numberOfCells[1]-1);
    }
    else if(domain.coordinates == 'spherical'){
        domain.xmin = domain.xmin;
        domain.xmax = domain.xmax;
        discretization.dlon = 60*(domain.xmax-domain.xmin)/(discretization.numberOfCells[0]-1);
        discretization.dlat = 60*(domain.ymax-domain.ymin)/(discretization.numberOfCells[1]-1);
    }

    // Initial condition
    earthquake = [];
    if(data.initialSurface){
        initialSurface = {
            array: data.initialSurface.array,
            textureId : 1,
            texture : undefined
        };
    }
    else if(data.earthquake){
        earthquake = data.earthquake;

        // longitude should be in [0,360]
        for(let i = 0; i<earthquake.length; i++){
            earthquake[i].ce = earthquake[i].ce;
        }
    }
            

    pcolorDisplay = {
        width : data.displayWidth,
        height : data.displayHeight
    };

    pois = output.pois ? output.pois : {};

    colors = [[ 0.01568627,  0.11372549,  0.23137255,  0.        ],
    [ 0.03137255,  0.23137255,  0.4627451 ,  0.01      ],
    [ 0.09411765,  0.30196078,  0.61568627,  0.011     ],
    [ 0.23137255,  0.41568627,  0.8       ,  0.05      ],
    [ 0.15294118,  0.1254902 ,  0.89411765,  0.051     ],
    [ 0.44313725,  0.72156863,  0.97647059,  0.1       ],
    [ 0.        ,  0.41568627,  0.06666667,  0.101     ],
    [ 0.        ,  0.81568627,  0.        ,  0.25      ],
    [ 0.5372549 ,  0.50980392,  0.        ,  0.251     ],
    [ 0.99607843,  0.89803922,  0.07843137,  0.5       ],
    [ 0.51372549,  0.31372549,  0.        ,  0.501     ],
    [ 0.88235294,  0.50196078,  0.0627451 ,  0.75      ],
    [ 0.62352941,  0.0745098 ,  0.        ,  0.751     ],
    [ 0.97647059,  0.10196078,  0.        ,  1.        ],
    [ 1.        ,  1.        ,  1.        ,  1.001     ],
    [ 1.        ,  0.25098039,  0.76862745,  2.25      ]];
    
    let cmin = -0.2;
    let cmax = 0.2;
    
    colors = [[ 0.        ,  0.        ,  0.3 ,  0.0       *(cmax-cmin)+cmin],
    [ 0.        ,  0.        ,  0.48666667,  0.06666667*(cmax-cmin)+cmin],
    [ 0.        ,  0.        ,  0.67333333,  0.13333333*(cmax-cmin)+cmin],
    [ 0.        ,  0.        ,  0.86      ,  0.2       *(cmax-cmin)+cmin],
    [ 0.06666667,  0.06666667,  1.        ,  0.26666667*(cmax-cmin)+cmin],
    [ 0.33333333,  0.33333333,  1.        ,  0.33333333*(cmax-cmin)+cmin],
    [ 0.6       ,  0.6       ,  1.        ,  0.4       *(cmax-cmin)+cmin],
    [ 0.86666667,  0.86666667,  1.        ,  0.46666667*(cmax-cmin)+cmin],
    [ 1.        ,  0.86666667,  0.86666667,  0.53333333*(cmax-cmin)+cmin],
    [ 1.        ,  0.6       ,  0.6       ,  0.6       *(cmax-cmin)+cmin],
    [ 1.        ,  0.33333333,  0.33333333,  0.66666667*(cmax-cmin)+cmin],
    [ 1.        ,  0.06666667,  0.06666667,  0.73333333*(cmax-cmin)+cmin],
    [ 0.9       ,  0.        ,  0.        ,  0.8       *(cmax-cmin)+cmin],
    [ 0.76666667,  0.        ,  0.        ,  0.86666667*(cmax-cmin)+cmin],
    [ 0.63333333,  0.        ,  0.        ,  0.93333333*(cmax-cmin)+cmin],
    [ 0.5       ,  0.        ,  0.        ,  1.0       *(cmax-cmin)+cmin]];
    
    colors = [[ 0.        ,  0.        ,  1.        ,  0.        *(cmax-cmin)+cmin],
    [ 0.13333333,  0.13333333,  1.        ,  0.06666667*(cmax-cmin)+cmin],
    [ 0.26666667,  0.26666667,  1.        ,  0.13333333*(cmax-cmin)+cmin],
    [ 0.4       ,  0.4       ,  1.        ,  0.2       *(cmax-cmin)+cmin],
    [ 0.53333333,  0.53333333,  1.        ,  0.26666667*(cmax-cmin)+cmin],
    [ 0.66666667,  0.66666667,  1.        ,  0.33333333*(cmax-cmin)+cmin],
    [ 0.8       ,  0.8       ,  1.        ,  0.4       *(cmax-cmin)+cmin],
    [ 0.93333333,  0.93333333,  1.        ,  0.46666667*(cmax-cmin)+cmin],
    [ 1.        ,  0.93333333,  0.93333333,  0.53333333*(cmax-cmin)+cmin],
    [ 1.        ,  0.8       ,  0.8       ,  0.6       *(cmax-cmin)+cmin],
    [ 1.        ,  0.66666667,  0.66666667,  0.66666667*(cmax-cmin)+cmin],
    [ 1.        ,  0.53333333,  0.53333333,  0.73333333*(cmax-cmin)+cmin],
    [ 1.        ,  0.4       ,  0.4       ,  0.8       *(cmax-cmin)+cmin],
    [ 1.        ,  0.26666667,  0.26666667,  0.86666667*(cmax-cmin)+cmin],
    [ 1.        ,  0.13333333,  0.13333333,  0.93333333*(cmax-cmin)+cmin],
    [ 1.        ,  0.        ,  0.        ,  1.        *(cmax-cmin)+cmin]];
    
    
    // flatten the array
    colors = colors.reduce((a,b)=>{
        return a.concat(b);
    });

    /* 
        Start WebGL context

    */

    let canvas = document.createElement("canvas");
    canvas.width = pcolorDisplay.width;
    canvas.height = pcolorDisplay.height;

    try
    {
        gl = canvas.getContext("webgl2");
        isWebGL2 = !!gl;
        if (!isWebGL2) {
            gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
        }
        gl.viewportWidth = canvas.width;
        gl.viewportHeight = canvas.height;
    } catch (e) {
    }
    if (!gl) {
        alert("Could not initialise Webgl.");
    }

    gl.enable(gl.DEPTH_TEST);

    gl.clearColor(0.0,0.0,0.0,1.0); //default color for any fbo, canvas included

    /*
        WebGL tools

    */

    let compileShader = function(type, source){
        let shader = gl.createShader(type); //shader handle
        gl.shaderSource(shader, source);
        gl.compileShader(shader);

        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            console.warn(source);
            throw gl.getShaderInfoLog(shader);
        }
        
        return shader;
    }

    let shaderProgram = function(vertexShader, fragmentShader){
        
        let uniforms = {}; // to store uniforms handles
        let program = gl.createProgram(); //program handle

        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        gl.linkProgram(program);
        
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            throw gl.getProgramInfoLog(program);
        }


        let uniformsCount = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);

        for(let i = 0; i< uniformsCount; i++){
            let uniformName = gl.getActiveUniform(program, i).name;
            if(uniformName.includes('[0]')) uniformName = uniformName.replace('[0]','');
            uniforms[uniformName] = gl.getUniformLocation(program, uniformName);
        }
        
        uniforms.vertexPositionAttribute = gl.getAttribLocation(program, "inPosition");
            
        
        gl.enableVertexAttribArray(uniforms.vertexPositionAttribute);
        
        return {uniforms,program};

    }

    let createShaders = function(){
        vertexShader = compileShader(gl.VERTEX_SHADER,`
            precision highp  float;
            attribute vec2 inPosition;

            varying  vec2 vUv;

            void main()
            {
                vUv = inPosition.xy*0.5+0.5;
                
                gl_Position = vec4(inPosition,0, 1);
            }    
        `);

        initialShader = compileShader(gl.FRAGMENT_SHADER,`
            precision highp float;
            
            uniform sampler2D bathymetry;
            uniform sampler2D initialSurface;

            varying vec2 vUv;

            void main()
            { 
                float h = texture2D(bathymetry, vUv).r;
                h = max(h,0.0);
                float eta = texture2D(initialSurface, vUv).r;
                gl_FragColor  = vec4(eta, 0.0, 0.0, h);
            }    
        `);
        
        okadaShader = compileShader(gl.FRAGMENT_SHADER,`
            precision highp float;
            varying vec2 vUv;
            uniform sampler2D bathymetry;
            
            uniform vec2 texel; 
            
            uniform float xmin;
            uniform float xmax;
            uniform float ymin;
            uniform float ymax;
            
            uniform sampler2D previousTexture;
            uniform float L;
            uniform float W;
            uniform float depth;
            uniform float slip;
            uniform float strike;
            uniform float dip;
            uniform float rake;
            uniform float U3;
            uniform float cn;
            uniform float ce;
            uniform int reference;
            uniform int coordinates;
            
            const float MY_PI = 3.14159265358979323;
            const float Rearth = 6378000.0;
            // const float Rearth = 6384.e+3;
            const int CENTER = 0;
            const int BEGINTOP = 1;
            const int MIDTOP = 2;
            const int BEGINBOTTOM = 3;
            const int MIDBOTTOM = 4;

            const int CARTESIAN = 0;
            const int SPHERICAL = 1;

            float I4(float db, float eta, float q, float dip, float nu, float R){
                float I = 0.0;
                if (abs(cos(dip))>1e-24){
                    I = (1.0-2.0*nu)/cos(dip)*(log(R+db)-sin(dip)*log(R+eta));
                }
                else{
                    I = -(1.0-2.0*nu)*q/(R+db);
                }
                return I;

            }

            float I5(float xi, float eta, float q, float dip,
                float nu, float R, float db){
                float x = sqrt(xi*xi+q*q);
                float I = 0.0;
                if (abs(cos(dip))>1e-24){
                    if(abs(xi)<1e-10){
                        I = 0.0;
                    }
                    else{
                        I = (1.0-2.0*nu)*2.0/cos(dip)*
                            atan( 
                                (  eta*(x+q*cos(dip)) + x*(R+x)*sin(dip)  )
                                /(  xi*(R+x)*cos(dip)   )
                            );                        
                    }
                
                }
                else{
                    I = -(1.0-2.0*nu)*xi*sin(dip)/(R+db);
                }
                return I;
            }

            float uz_ss(float xi, float eta, float q, float dip, float nu){
                float R = sqrt(xi*xi+eta*eta+q*q);
                float db = eta*sin(dip) - q*cos(dip);
                float u = db*q/(R*(R+eta))
                        + q*sin(dip)/(R+eta)
                        + I4(db,eta,q,dip,nu,R)*sin(dip);
                return u;
            }

            float uz_ds(float xi, float eta, float q, float dip, float nu){
                float R = sqrt(xi*xi+ eta*eta+ q*q);
                float db = eta*sin(dip) - q*cos(dip);
                float u = db*q/(R*(R+xi))
                    - I5(xi,eta,q,dip,nu,R,db)*sin(dip)*cos(dip);
                if (q*R!=0.0){
                    u = u + sin(dip)*atan(xi*eta/(q*R));
                }
                else{
                    if(xi*eta>0.0){
                        u = u + MY_PI*sin(dip);
                    }
                    else if(xi*eta<0.0){
                        u = u - MY_PI*sin(dip);
                    }
                }
                return u;
            }

            float uz_tf(float xi, float eta, float q, float dip, float nu){
                float R = sqrt(xi*xi+eta*eta+q*q);
                float db = eta*sin(dip) - q*cos(dip);
                float u = (eta*cos(dip)+q*sin(dip))*q/(R*(R+xi))
                        + cos(dip)*xi*q/(R*(R+eta))
                        - I5(xi,eta,q,dip,nu,R,db)*sin(dip)*sin(dip);
                if (q*R!=0.0){
                    u = u - cos(dip)*atan(xi*eta/(q*R));
                }
                else{
                    if(xi*eta>0.0){
                        u = u - MY_PI*cos(dip);
                    }
                    else if(xi*eta<0.0){
                        u = u + MY_PI*cos(dip);
                    }
                }
                return u;
            }

            //finally
            float okada(float n, float e, float depth, float strike, float dip, float L,
                        float W, float rake, float slip, float U3){
                float nu = 0.25;
                float pi = 3.14159265359;
                strike = strike*pi/180.0;
                dip = dip*pi/180.0;
                rake = rake*pi/180.0;

                float U1 = cos(rake)*slip;
                float U2 = sin(rake)*slip;
                
                // const int CENTER = 0;
                // const int BEGINTOP = 1;
                // const int MIDTOP = 2;
                // const int BEGINBOTTOM = 3;
                // const int MIDBOTTOM = 4;
        
                // rotate (n,e) coordinates centered at the reference point (after stereo proj)
                float x = cos(strike)*n + sin(strike)*e;
                float y = sin(strike)*n - cos(strike)*e;
                float d = depth;

                // translate according to reference point location
                // begin bottom is the default case for this formula
                if(reference==CENTER){
                    x = x + 0.5*L;
                    y = y + 0.5*W*cos(dip);
                    d = d + 0.5*W*sin(dip);         
                }
                else if(reference == BEGINTOP){
                    y = y + W*cos(dip);
                    d = d + W*sin(dip);
                }
                else if(reference == MIDTOP){
                    x = x + 0.5*L;
                    y = y + W*cos(dip);
                    d = d + W*sin(dip);
                }
                else if(reference == MIDBOTTOM){
                    x = x + 0.5*L;
                }


                float p = y*cos(dip) + d*sin(dip);
                float q = y*sin(dip) - d*cos(dip);

                float F1 = uz_ss(x  ,p  ,q,dip,nu);
                float F2 = uz_ss(x  ,p-W,q,dip,nu);
                float F3 = uz_ss(x-L,p  ,q,dip,nu);
                float F4 = uz_ss(x-L,p-W,q,dip,nu);
                float F = F1-F2-F3+F4;

                float G1 = uz_ds(x  ,p  ,q,dip,nu);
                float G2 = uz_ds(x  ,p-W,q,dip,nu);
                float G3 = uz_ds(x-L,p  ,q,dip,nu);
                float G4 = uz_ds(x-L,p-W,q,dip,nu);
                float G = G1-G2-G3+G4;

                float H1 = uz_tf(x  ,p  ,q,dip,nu);
                float H2 = uz_tf(x  ,p-W,q,dip,nu);
                float H3 = uz_tf(x-L,p  ,q,dip,nu);
                float H4 = uz_tf(x-L,p-W,q,dip,nu);
                float H = H1-H2-H3+H4;


                float uz = -U1/(2.0*pi)*F - U2/(2.0*pi)*G + U3/(2.0*pi)*H ;

                return uz;
            }

            // Stereographic projection function

            vec2 stereographic_projection(float latin, float lonin,
                                            float lat0, float lon0){
                /*
                    Gives the stereographic projection of points (lat,lon)
                    onto the plane tangent to the Earth's ellipsoid in (lat0,lon0)
                    Input:
                        latin,lonin : coordinates of points in degrees (array)
                        lat0,lon0: coordinates of tangency point in degrees (float)
                */
                float pi = 3.141592653589793 ;
                float pole = pi/2.0 - 1e-5;
                float rad_deg = pi/180.0;

                float lat = latin*rad_deg;
                float lon = lonin*rad_deg;
                float lt0 = lat0*rad_deg;
                float ln0 = lon0*rad_deg;

                lat = max(-pole,min(pole,lat));
                lt0 = max(-pole,min(pole,lt0));
                float CS = cos(lat);
                float SN = sin(lat);
                float CS0 = cos(lt0);
                float SN0 = sin(lt0);
                float xf = 0.0; //false easting
                float yf = 0.0; //false northing

                float A = 6378137.0000;            //ELLIPSOIDAL SEMI-MAJOR AXIS
                float B = 6356752.3142;            //ELLIPSOIDAL SEMI-MINOR AXIS
                float F = 0.00335281067183;      // FLATTENING, F = (A-B)/A
                float E = 0.08181919092891;        // ECCENTRICITY, E = SQRT(2.0*F-F**2)
                float F2 = 0.00669438000426;        // F2 = E**2
                float ES = 0.00673949675659;        // 2ND ECCENTRICITY, ES = E**2/(1-E**2)

                float K0 = 0.9996;                // SCALE FACTOR

                float TMP = sqrt(1.0-F2*SN*SN);
                float TMP0 = sqrt(1.0-F2*SN0*SN0);
                float RHO0 = A*(1.0-F2)/(TMP0*TMP0*TMP0);
                float NU0 = A/TMP0;
                float R = sqrt(RHO0*NU0);
                float N = sqrt(1.0+F2*CS0*CS0*CS0*CS0/(1.0-F2));

                float S1 = (1.0+SN0)/(1.0-SN0);
                float S2 = (1.0-E*SN0)/(1.0+E*SN0);
                float W1 = pow((S1*pow(S2,E)),N);
                float SN_XI0 = (W1-1.0)/(W1+1.0);
                float C = (N+SN0)*(1.0-SN_XI0)/(N-SN0)/(1.0+SN_XI0);

                float W2 = C*W1;
                float SA = (1.0+SN)/(1.0-SN);
                float SB = (1.0-E*SN)/(1.0+E*SN);
                float W = pow(C*(SA*pow(SB,E)),N);


                float XI0 = asin((W2-1.0)/(W2+1.0));
                float LM0 = ln0;

                float LM = N*(lon-LM0)+LM0;
                float XI = asin((W-1.0)/(W+1.0));

                float BETA = 1.0 + sin(XI)* sin(XI0) + cos(XI)*cos(XI0)*cos(LM-LM0);

                float Y = yf + 2.0*R*K0*(sin(XI)*cos(XI0)
                                - cos(XI)*sin(XI0)*cos(LM-LM0))/BETA;
                float X = xf + 2.0*R*K0*cos(XI)*sin(LM-LM0)/BETA;

                vec2 XY = vec2(X,Y);
                return XY;
            }

            vec2 simpleProjection(float latin, float lonin, float lat0, float lon0){
                float pi = 3.141592653589793 ;
                float y = Rearth*(latin-lat0)*pi/180.0;
                float x = Rearth*cos(lat0*pi/180.0)*(lonin-lon0)*pi/180.0;

                return vec2(x,y);

            }

            void main()
            {

                // normalize vUv so it is defined 1-1 in [0,1] 
                // this way the first pixel corresponds to xmin and the last to xmax, exactly

                float nx = 1.0/texel.x;
                float ny = 1.0/texel.y;
                float V = (vUv.y-0.5*texel.y)/((ny-1.0)*texel.y);
                float U = (vUv.x-0.5*texel.x)/((nx-1.0)*texel.x);
                float n = ymin + V*(ymax-ymin);
                float e = xmin + U*(xmax-xmin);

                // center on reference point and make projection if necessary
                vec2 pos;
                if(coordinates==CARTESIAN){
                    pos = vec2(e-ce,n-cn);
                }
                else if(coordinates==SPHERICAL){
                    // vec2 pos = stereographic_projection(n,e,cn,ce);
                    pos = simpleProjection(n,e,cn,ce);
                }
                
                // only calculate in a square around the rectangle
                float value = 0.0;
                if (abs(pos.x)<L/2.0*8.0 && abs(pos.y)<L/2.0*8.0){

                    value  = okada(pos.g, pos.r, depth, strike, dip, L, W, rake, slip, U3);

                }
                
                float bathymetry = texture2D(bathymetry, vUv).r;
                value = value*step(0.0,bathymetry);
                bathymetry = max(0.0, bathymetry);

                float u0 = texture2D(previousTexture, vUv).r;

                value = value + u0;


                gl_FragColor = vec4(value,0.0,0.0, bathymetry+value);

            }
        
        `);
        cartesianWaveShader = compileShader(gl.FRAGMENT_SHADER,`
            precision highp float;
            
            uniform sampler2D u0;
            
            uniform vec2 texel;

            varying vec2 vUv;

            const float g = 9.81;
            const float eps = 1e-2;

            uniform float dt;
            uniform float dx;
            uniform float dy;
            

            float massEquation(vec4 uij, vec4 uimj, vec4 uijm){            

                float eta = 0.0;
                float hij = uij.a;
                if(hij > eps){
                    eta = uij.r -dt/dx*(uij.g-uimj.g) -dt/dy*(uij.b-uijm.b);
                }

                return  eta;

            }

            void main()
            { 
                vec4 uij  = texture2D(u0, vUv );
                vec4 uipj = texture2D(u0, vUv+vec2(texel.x,0.0));
                vec4 uijp = texture2D(u0, vUv+vec2(0.0,texel.y));
                vec4 uimj = texture2D(u0, vUv+vec2(-texel.x,0.0));
                vec4 uijm = texture2D(u0, vUv+vec2(0.0,-texel.y));
                
                vec4 uipjm = texture2D(u0, vUv+vec2(texel.x,-texel.y));
                vec4 uimjp = texture2D(u0, vUv+vec2(-texel.x,texel.y));

                vec4 u2ij = vec4(0.0);
                vec4 u2ipj = vec4(0.0);
                vec4 u2ijp = vec4(0.0);

                u2ij.r = massEquation(uij, uimj, uijm);

                u2ipj.r = massEquation(uipj, uij, uipjm);

                u2ijp.r = massEquation(uijp, uimjp, uij);

                float hij = uij.a;
                float hipj = uipj.a;
                float hijp = uijp.a;

                if(hij>eps && hipj>eps){
                    // if not dry, otherwise defaults to 0.0;
                    u2ij.g = uij.g - g*hij*dt/dx*(u2ipj.r - u2ij.r);
                }
                
                if(hij>eps && hijp>eps){
                    u2ij.b = uij.b - g*hij*dt/dy*(u2ijp.r - u2ij.r);
                }

                u2ij.a = uij.a;
                
                gl_FragColor  = u2ij;
            }    
        `);

        sphericalWaveShader = compileShader(gl.FRAGMENT_SHADER,`
            precision highp float;

            varying vec2 vUv;
            uniform sampler2D u0; 
            uniform vec2 texel; 

            uniform float dlon;
            uniform float dlat;
            uniform float dt;
            uniform float xmin; 
            uniform float xmax;
            uniform float ymin;
            uniform float ymax;
            uniform int isPeriodic;

            const float rad_min = 0.000290888208665721; 
            const float rad_deg = 0.01745329252;
            const float cori_w = 7.2722e-5;
            const float gx = 1e-5;
            const float g = 9.81;
            const float Rearth = 6378000.0;


            float degToRad(float degrees){
                // convert from degrees to radians
                return rad_deg*degrees;
            }

            float minToRad(float minutes){
                // convert from minutes to radians
                return rad_min*minutes;
            }

            float openBoundary(vec2 vUv, vec4 u_ij, vec4 u_ijm, vec4 u_imj, float h_ij){
                float eta;
            
                float etaij = u_ij.r;
                float Mij = u_ij.g;
                float Nij = u_ij.b;
                
                float etaimj = u_imj.r;
                float Mimj = u_imj.g;
                float Nimj = u_imj.b;
            
                float etaijm = u_ijm.r;
                float Mijm = u_ijm.g;
                float Nijm = u_ijm.b;
                
                //j=0
                float k = 1.0;
                if (vUv.y <= k*texel.y){
                    eta = 0.0;
                    if (h_ij>gx){
                        float c = sqrt(g*h_ij);
                        float z = sqrt(Nij*Nij+0.25*(Mij+Mimj)*(Mij+Mimj))/c;
                        if (Nij>0.0){
                            z = -z;
                        }
                        eta = z;
            
                    }		
                }
            
                if (vUv.y >= 1.0-k*texel.y){
                    eta = 0.0;
                    if (h_ij>gx){
                        float c = sqrt(g*h_ij);
                        float z = sqrt(Nijm*Nijm+0.25*(Mij+Mimj)*(Mij+Mimj))/c;
                        if (Nijm<0.0){
                            z = -z;
                        }
                        eta = z;
                    }
                }
            
                if(xmax-xmin<360.0-0.2){
                    if (vUv.x <= k*texel.x){
                        eta = 0.0;
                        if (h_ij>gx){
                            float c = sqrt(g*h_ij);
                            float z = sqrt(Mij*Mij+0.25*(Nij+Nijm)*(Nij+Nijm))/c;
                            if (Mij>0.0){
                                z = -z;
                            }
                            eta = z;
            
                        }
                        
                    }
                    
                    if (vUv.x >=1.0-k*texel.x){
                        eta = 0.0;
                        if (h_ij>gx){
                            float c = sqrt(g*h_ij);
                            float z = sqrt(Mimj*Mimj+0.25*(Nij+Nijm)*(Nij+Nijm))/c;
                            if (Mimj<0.0){
                                z = -z;
                            }
                            eta = z;
            
                        }
                        
                    }
                }
                
            
                if(vUv.x <= k*texel.x && vUv.y<=k*texel.y){
                    eta = 0.0;
                    if (h_ij>gx){
                        float c = sqrt(g*h_ij);
                        float z = sqrt(Mij*Mij +Nij*Nij)/c;
                        if (Nij>0.0){
                            z = -z;
                        }
                        eta = z;
                    }
                }
            
                if(vUv.x >= 1.0-k*texel.x && vUv.y<=k*texel.y){
                    eta = 0.0;
                    if (h_ij>gx){
                        float c = sqrt(g*h_ij);
                        float z = sqrt(Mimj*Mimj+Nij*Nij)/c;
                        if (Nij>0.0){
                            z = -z;
                        }
                        eta = z;
                    }
                }
            
                if(vUv.x <= k*texel.x && vUv.y>=1.0-k*texel.y){
                    eta = 0.0;
                    if (h_ij>gx){
                        float c = sqrt(g*h_ij);
                        float z = sqrt(Mij*Mij +Nij*Nij)/c;
                        if (Nijm<0.0){
                            z = -z;
                        }
                        eta = z;
                    }
                }
            
                if(vUv.x >= 1.0 - k*texel.x && vUv.y>=1.0-k*texel.y){
                    eta = 0.0;
                    if (h_ij>gx){
                        float c = sqrt(g*h_ij);
                        float z = sqrt(Mimj*Mimj +Nijm*Nijm)/c;
                        if (Nijm<0.0){
                            z = -z;
                        }
                        eta = z;
                    }
                }
                
                return eta;
            }            

            float updateInnerCellSurface(vec2 vUv, vec2 UV, vec4 uij, vec4 uimj, vec4 uijm){
                /* apply mass conservation equation to update water surface */

                // calculate space varying constants
                float V  = UV.y;
                float dlonrad = minToRad(dlon);
                float dlatrad = minToRad(dlat);
                float latj = V*(ymax-ymin)+ymin;
                float latjp = latj + 0.5*dlat/60.0;
                float latjm = latj - 0.5*dlat/60.0;
                float coslatj = cos(degToRad(latj));
                float coslatjp = cos(degToRad(latjp));
                float coslatjm = cos(degToRad(latjm));

                float eta2ij;
            
                eta2ij = uij.r - dt/(Rearth*coslatj*dlonrad)*(uij.g - uimj.g + uij.b*coslatjp - uijm.b*coslatjm); 


                return eta2ij;

            }

            float updateSurface(vec2 vUv, vec2 UV, vec4 uij, vec4 uimj, vec4 uijm){

                bool isBoundary = vUv.x < texel.x || vUv.x > 1.0 - texel.x || vUv.y < texel.y || vUv.y > 1.0 - texel.y;

                float eta2ij = 0.0;
                
                if (isPeriodic == 0 && isBoundary){
                    
                    eta2ij = openBoundary(UV, uij, uijm, uimj, uij.a);


                }
                else if (uij.a > gx) {

                    eta2ij = updateInnerCellSurface(vUv, UV, uij, uimj, uijm);

                }               

                return eta2ij;

            }
            vec2 correctedUV(vec2 uv){
                // corrected uv coordinates 
                // so min(u) maps to U=0, and max(u) maps to U = 1
                // same with v
                // since normally, these are displace by half texel

                float nx = 1.0/texel.x;
                float ny = 1.0/texel.y;
                float V = (uv.y-0.5*texel.y)/((ny-1.0)*texel.y);
                float U = (uv.x-0.5*texel.x)/((nx-1.0)*texel.x);

                return vec2(U,V);
            }
            void main()
            {
                // u = (eta, p, q, h)
                // eta: free surface
                // p: x-momentum
                // q: y-momentum
                // h: water depth, >0 if wet, <0 if dry.
                

                // read previous frame, handling periodic boundaries
                vec2 right = vec2(texel.x,0.0);
                vec2 front = vec2(0.0, texel.y);
                vec4 uij = texture2D(u0, vUv);                
                
                vec4 uimj, uipj, uimjp, uipjm;

                if(isPeriodic == 1){
                    float uright = mod(vUv.x + right.x, 1.0);
                    float uleft = mod(vUv.x - right.x, 1.0);
                    
                    uimj = texture2D(u0, vec2(uleft, vUv.y));
                    uipj = texture2D(u0, vec2(uright, vUv.y));
                    uimjp = texture2D(u0, vec2(uleft, vUv.y) + front);
                    uipjm = texture2D(u0, vec2(uright, vUv.y) - front);
                }
                else{
                    uimj = texture2D(u0, vUv - right);
                    uipj = texture2D(u0, vUv + right);
                    uimjp = texture2D(u0, vUv - right + front);
                    uipjm = texture2D(u0, vUv + right - front);

                }

                vec4 uijp = texture2D(u0, vUv + front);
                vec4 uijm = texture2D(u0, vUv - front);
                
                // mass conservation
                float eta2ij = updateSurface(vUv, correctedUV(vUv), uij, uimj, uijm);

                float eta2ipj = updateSurface(vUv+right, correctedUV(vUv+right), uipj, uij, uipjm);

                float eta2ijp = updateSurface(vUv+front, correctedUV(vUv+front), uijp, uimjp, uij);

                // momentum conservation 

                float hiPlusHalfj = 0.5*(uij.a + uipj.a);
                float hijPlusHalf = 0.5*(uij.a + uijp.a);
                vec2 UV = correctedUV(vUv);
                float U = UV.x;
                float V  = UV.y;
                float lon = U*(ymax-ymin)+ymin;
                float latj = V*(ymax-ymin)+ymin;
                float coslatj = cos(degToRad(latj));
                

                float M2ij = 0.0;
                
                if(hiPlusHalfj > gx){
                    M2ij = uij.g - dt*g*hiPlusHalfj/(Rearth*coslatj*minToRad(dlon))*(eta2ipj - eta2ij);
                }

                float N2ij = 0.0;
                
                if(hijPlusHalf > gx){
                    N2ij = uij.b - dt*g*hijPlusHalf/(Rearth*minToRad(dlat))*(eta2ijp - eta2ij);
                }

            
                gl_FragColor = vec4( eta2ij, M2ij, N2ij, uij.a);

            }

        `);

        maxHeightsShader = compileShader(gl.FRAGMENT_SHADER,`
            precision highp float;

            uniform sampler2D maxHeights;
            uniform sampler2D wave;
            uniform float currentTime;

            varying vec2 vUv;

            void main(){
                vec4 pixelMaxHeight = texture2D(maxHeights, vUv);
                float maxHeight = pixelMaxHeight.r;
                float arrivalTime = pixelMaxHeight.g;

                float newHeight = texture2D(wave, vUv).r;

                if( maxHeight < newHeight){
                    maxHeight = newHeight;
                }
                
                if(arrivalTime == 0.0 && newHeight > 1e-5){
                    arrivalTime = currentTime;
                }
                gl_FragColor = vec4( maxHeight, arrivalTime, 0.0, 1.0);

            }
        `);
        
        displayShader = compileShader(gl.FRAGMENT_SHADER,`
            precision highp float;
            
            const int nColors = 16;            
            uniform sampler2D field;
            uniform vec4 colormap[16];
            uniform int displayedChannel;
            
            varying vec2 vUv;

            vec3 getPseudoColor(float value){
                vec3 pseudoColor;

                if(value <= colormap[0].a){
                    pseudoColor = colormap[0].rgb;
                }
                else if (value > colormap[16-1].a){
                    pseudoColor = colormap[16-1].rgb;
                }
                else{
                    for (int i=1; i<16; i++){
                        vec4 cleft = colormap[i-1];
                        vec4 cright = colormap[i];
            
                        if (value>cleft.a && value <=cright.a){
                            float t = (value - cleft.a)/(cright.a - cleft.a);
                            pseudoColor = mix(cleft.rgb, cright.rgb, t);
                            break;
                        }
                    }
                }
            
                return pseudoColor;
            }
            
            void main()
            { 
                float uij  = texture2D(field, vUv).r;
                if(displayedChannel == 1){
                    uij = texture2D(field, vUv).g;
                }
                else if(displayedChannel == 2){
                    uij = texture2D(field, vUv).b;
                }
                else if(displayedChannel == 3){
                    uij = texture2D(field, vUv).a;
                }

                vec3 color = getPseudoColor(uij);

                float alpha  = pow(abs(uij),0.2);
                gl_FragColor  = vec4(color, alpha);
            }    
        `);        
        initialProgram = shaderProgram(vertexShader, initialShader);
        okadaProgram = shaderProgram(vertexShader, okadaShader);
        cartesianWaveProgram = shaderProgram(vertexShader, cartesianWaveShader);
        sphericalWaveProgram = shaderProgram(vertexShader, sphericalWaveShader);
        displayProgram = shaderProgram(vertexShader, displayShader);
        maxHeightsProgram = shaderProgram(vertexShader, maxHeightsShader);
    }

    let createBuffers = function(){ 
        let vertexPositionBufferHandle = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, vertexPositionBufferHandle);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, -1, 1, 1, 1, 1, -1]), gl.STATIC_DRAW);

        let facesBufferHandle = gl.createBuffer();
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, facesBufferHandle);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array([0, 1, 2, 0, 2, 3]),gl.STATIC_DRAW);
        
        gl.vertexAttribPointer(cartesianWaveProgram.uniforms.vertexPositionAttribute, 2, gl.FLOAT, false, 0, 0);
    }

    let createTextureFromData = function( width, height, data, textureId, internalFormat, format, type ) {
        // creates a texture from an array "data" of size "width"*height"*4 with the given webgl formats and id

        if (!gl.getExtension("OES_texture_float")&& !gl.getExtension("EXT_color_buffer_float")) {
            throw("Requires OES_texture_float extension");
        }
        
        var texture = gl.createTexture();
        gl.activeTexture(gl.TEXTURE0 + textureId);
        gl.bindTexture( gl.TEXTURE_2D, texture );

        // 32F: diferencia con webgl/webgl2
        
        gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, width, height, 0, format, type, new Float32Array(data) );

        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE );
        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE );
        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST );
        gl.texParameteri( gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST );
        return {texture, textureId};
            
    }

    let createTextureFromMatrix = function(matrix,textureId){

        let internalFormat = isWebGL2? gl.RGBA32F : gl.RGBA;   
        let format = gl.RGBA;
        let type = gl.FLOAT;

        let raveledMatrix = [];
        for(let j =0; j<matrix.length;j++){
            for(let i=0;i<matrix[0].length; i++){
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
            textureId, internalFormat, format, type )

        return texture;    
    }

    let createFBO = function(textureId, w, h,internalFormat, format, type, param ){
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

        gl.activeTexture(gl.TEXTURE0+textureId);
        let texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, param);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, param);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, w, h, 0, format, type, null); //data will be rendered later

        let fbo = gl.createFramebuffer();
        gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
        gl.viewport(0, 0, w, h);
        gl.clear(gl.COLOR_BUFFER_BIT);
        return {texture, fbo, textureId};
    }

    let createDoubleFBO = function(textureId1, textureId2, w, h,internalFormat, format, type, param ){
        
        let fbo1 = createFBO(textureId1, w, h,internalFormat, format, type, param )
        let fbo2 = createFBO(textureId2, w, h,internalFormat, format, type, param )

        return {
            get first(){
                return fbo1
            },
            get second(){
                return fbo2
            },
            swap: function(){
                let temp = fbo1;
                fbo1 = fbo2;
                fbo2 = temp;
            }
        }
    }


    let renderFrameBuffer = function(frameBuffer){
        gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuffer);
        gl.drawElements(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0);
    }

    let renderInitialProgram = function(){
        
        gl.viewport(0, 0, discretization.numberOfCells[0], discretization.numberOfCells[1]);
        gl.useProgram(initialProgram.program);
        gl.uniform1i(initialProgram.uniforms.bathymetry, bathymetry.texture.textureId);
        gl.uniform1i(initialProgram.uniforms.initialSurface, initialSurface.texture.textureId);

        renderFrameBuffer(wave.first.fbo);
    }

    let renderOkadaProgram = function(finiteFault){
        
        gl.viewport(0, 0, discretization.numberOfCells[0], discretization.numberOfCells[1]);
        gl.useProgram(okadaProgram.program);
        gl.uniform1i(okadaProgram.uniforms.bathymetry, bathymetry.texture.textureId);


        gl.uniform2f(okadaProgram.uniforms.texel, 1/discretization.numberOfCells[0], 1/discretization.numberOfCells[1]);

        gl.uniform1f(okadaProgram.uniforms.xmin, domain.xmin) ;
        gl.uniform1f(okadaProgram.uniforms.xmax, domain.xmax) ;
        gl.uniform1f(okadaProgram.uniforms.ymin, domain.ymin) ;
        gl.uniform1f(okadaProgram.uniforms.ymax, domain.ymax) ;


        gl.uniform1f(okadaProgram.uniforms.L, finiteFault.L);
        gl.uniform1f(okadaProgram.uniforms.W, finiteFault.W);
        gl.uniform1f(okadaProgram.uniforms.depth, finiteFault.depth);
        gl.uniform1f(okadaProgram.uniforms.slip, finiteFault.slip);
        gl.uniform1f(okadaProgram.uniforms.strike, finiteFault.strike);
        gl.uniform1f(okadaProgram.uniforms.dip, finiteFault.dip);
        gl.uniform1f(okadaProgram.uniforms.rake, finiteFault.rake);
        gl.uniform1f(okadaProgram.uniforms.U3, finiteFault.U3);
        gl.uniform1f(okadaProgram.uniforms.cn, finiteFault.cn);
        gl.uniform1f(okadaProgram.uniforms.ce, finiteFault.ce);
        gl.uniform1i(okadaProgram.uniforms.previousTexture, wave.first.textureId);


        if(finiteFault.reference == 'center'){
            gl.uniform1i(okadaProgram.uniforms.reference, 0);
        }
        else if(finiteFault.reference == 'begin top'){
            gl.uniform1i(okadaProgram.uniforms.reference, 1);
        }
        else if(finiteFault.reference == 'mid top'){
            gl.uniform1i(okadaProgram.uniforms.reference, 2);
        }
        else if(finiteFault.reference == 'begin bottom'){
            gl.uniform1i(okadaProgram.uniforms.reference, 3);
        }
        else if(finiteFault.reference == 'mid bottom'){
            gl.uniform1i(okadaProgram.uniforms.reference, 4);
        }
        
        if(domain.coordinates == 'cartesian'){
            gl.uniform1i(okadaProgram.uniforms.coordinates, 0);
        }
        else if(domain.coordinates == 'spherical'){
            gl.uniform1i(okadaProgram.uniforms.coordinates, 1);
        }

        renderFrameBuffer(wave.second.fbo);
    }

    let renderEarthquake = () => {
        for(let i = 0; i < earthquake.length; i++){
            renderOkadaProgram(earthquake[i]);
            wave.swap();
        }
    }

    let renderCartesianProgram = function(){

        gl.viewport(0, 0, discretization.numberOfCells[0], discretization.numberOfCells[1]);
        gl.useProgram(cartesianWaveProgram.program);
        gl.uniform2f(cartesianWaveProgram.uniforms.texel, 1/discretization.numberOfCells[1], 1/discretization.numberOfCells[1]);
        gl.uniform1i(cartesianWaveProgram.uniforms.u0, wave.first.textureId);

        gl.uniform1f(cartesianWaveProgram.uniforms.dx, discretization.dx);
        gl.uniform1f(cartesianWaveProgram.uniforms.dy, discretization.dy);
        gl.uniform1f(cartesianWaveProgram.uniforms.xmin, domain.xmin);
        gl.uniform1f(cartesianWaveProgram.uniforms.xmax, domain.xmax);
        gl.uniform1f(cartesianWaveProgram.uniforms.ymin, domain.ymin);
        gl.uniform1f(cartesianWaveProgram.uniforms.ymax, domain.ymax);
        gl.uniform1f(cartesianWaveProgram.uniforms.dt, discretization.dt);

        renderFrameBuffer(wave.second.fbo);
        wave.swap();

    }

    let renderSphericalProgram = function(){        
        gl.useProgram(sphericalWaveProgram.program);
        gl.viewport(0, 0, discretization.numberOfCells[0], discretization.numberOfCells[1]);


        gl.uniform1f(sphericalWaveProgram.uniforms.dlon, discretization.dlon);
        gl.uniform1f(sphericalWaveProgram.uniforms.dlat, discretization.dlat);
        gl.uniform1f(sphericalWaveProgram.uniforms.dt, discretization.dt);
        gl.uniform1f(sphericalWaveProgram.uniforms.xmin, domain.xmin);
        gl.uniform1f(sphericalWaveProgram.uniforms.xmax, domain.xmax);
        gl.uniform1f(sphericalWaveProgram.uniforms.ymin, domain.ymin);
        gl.uniform1f(sphericalWaveProgram.uniforms.ymax, domain.ymax);

        
        gl.uniform2f(sphericalWaveProgram.uniforms.texel, 1/discretization.numberOfCells[0], 1/discretization.numberOfCells[1]);
        gl.uniform1i(sphericalWaveProgram.uniforms.u0, wave.first.textureId);
        gl.uniform1i(sphericalWaveProgram.uniforms.isPeriodic, domain.isPeriodic);
        renderFrameBuffer(wave.second.fbo);
        wave.swap();

    }

    let renderMaxHeightsProgram = () => {
        gl.useProgram(maxHeightsProgram.program);
        gl.viewport(0, 0, discretization.numberOfCells[0], discretization.numberOfCells[1]);

        gl.uniform1f(maxHeightsProgram.uniforms.currentTime, discretization.stepNumber * discretization.dt);
        gl.uniform1i(maxHeightsProgram.uniforms.maxHeights, maxHeights.first.textureId);
        gl.uniform1i(maxHeightsProgram.uniforms.wave, wave.first.textureId);
        renderFrameBuffer(maxHeights.second.fbo);
        maxHeights.swap();
    }
        
    let renderDisplayProgram = function(){
        
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        gl.useProgram(displayProgram.program);
        gl.uniform4fv(displayProgram.uniforms.colormap, new Float32Array(colors));
        
        let displayedChannel = 0 ;
        if( displayOption === 'heights'){
            gl.uniform1i(displayProgram.uniforms.field, wave.first.textureId );//TDDO: fix texid
        }
        else if( displayOption === 'max heights' || displayOption === 'arrival times'){

            gl.uniform1i(displayProgram.uniforms.field, maxHeights.first.textureId );
            if(displayOption === 'arrival times') displayedChannel = 1;
        }
        gl.uniform1i(displayProgram.uniforms.displayedChannel, displayedChannel)

        renderFrameBuffer(null);

    }

    let initFBOs = function(){
        const internalFormat = isWebGL2 ? gl.RGBA32F : gl.RGBA;
        const format = gl.RGBA; 
        const textype = gl.FLOAT;
        let textureIdHeights1 = 2;
        let textureIdHeights2 = 3;
        let textureIdMaxHeights1 = 4;
        let textureIdMaxHeights2 = 5;
        let param = gl.NEAREST;

        wave = createDoubleFBO(textureIdHeights1, textureIdHeights2, discretization.numberOfCells[0], discretization.numberOfCells[1], 
            internalFormat, format, textype, param);
        maxHeights = createDoubleFBO(textureIdMaxHeights1, textureIdMaxHeights2, discretization.numberOfCells[0], discretization.numberOfCells[1], 
            internalFormat, format, textype, param);

        if(initialSurface){
            renderInitialProgram();
        }
        else if(earthquake){
            renderEarthquake();
        }

        renderMaxHeightsProgram();

        renderDisplayProgram();
    }

    let readFBOPixels = function (frameBufferObject, left, top, width, height) {

        var pixelData = new Float32Array(width * height * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, frameBufferObject);
        gl.readPixels(left, top, width, height, gl.RGBA, gl.FLOAT, pixelData);
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        
        return pixelData;
    }

    let setPOIs = function(){
        /*
            Initializes the list of pois in the model
        */
        Object.keys(pois).forEach(function(poi){
            var dlon = discretization.dlon;
            var dlat = discretization.dlat;
            var lowerLeftCorner = [domain.xmin, domain.ymin];

            pois[poi].location[0] = pois[poi].location[0];

            var i = Math.floor((pois[poi].location[0]-lowerLeftCorner[0])/(dlon/60.0)+0.5);
            var j = Math.floor((pois[poi].location[1]-lowerLeftCorner[1])/(dlat/60.0)+0.5);

            pois[poi].pixel = [i+1,j+1];
            pois[poi].surface = [];
            pois[poi].time = [];
        });
    }

    let storePOISValues = function(){
        
        // save current time in minutes

        // poisTime.push(simulationData.currentIterationTime);

        // store water surface elevation at points of interest (POIs)

        Object.keys(pois).forEach(function(poi){
            var pixelsurface = readFBOPixels(wave.first.fbo, pois[poi].pixel[0], pois[poi].pixel[1], 1, 1)[0];

            pois[poi].surface.push(pixelsurface);
            pois[poi].time.push(discretization.dt*discretization.stepNumber);
            
        });
    }

    let exportBuffer = function(fbo, variableIndex=0, iStart=0, jStart = 0, 
         Lx = discretization.numberOfCells[0], Ly = discretization.numberOfCells[1]){
        let array = readFBOPixels(fbo, iStart, jStart, Lx, Ly);
        array = array.filter((elem,index)=>{
            return (index-variableIndex) % 4 == 0;
        });
        
        return array;
    }

    let nextgrid = 0;
    let runSimulationStep = function(){
        discretization.stepNumber++;

        if(domain.coordinates == 'cartesian'){
            renderCartesianProgram();
            
        }
        else if(domain.coordinates == 'spherical'){
            renderSphericalProgram();                    
        }

        renderMaxHeightsProgram();

        storePOISValues();
              
    }

    let setTimeStep = function(options){
        let hmax = Math.max.apply(null, 
            bathymetry.array.map(
                (row)=>{return Math.max.apply(null,row)}
            )
        );

        let cfl;
        if(options.timeStep !== undefined){
            discretization.dt = options.timeStep;
        }
        else{
            /* If dt is not given use given or predefined cfl */ 

            cfl = options.cfl === undefined ? 0.5 : options.cfl;
            if(domain.coordinates == 'cartesian'){
                discretization.dt = cfl * Math.min(discretization.dx,discretization.dy)/Math.sqrt(hmax*9.81);
            }
            else if(domain.coordinates == 'spherical'){
                // let Rearth = 6378000.0;
                let Rearth = 6378000;
                let radPerDeg = 0.01745329252;
                let radPerMin = 0.000290888208665721;
                var latMax = Math.max(Math.abs(domain.ymax),Math.abs(domain.ymin));//Math.max(Math.abs(ymin),Math.abs(ymax));
                var dxReal = Rearth * Math.cos(latMax * radPerDeg) * discretization.dlon * radPerMin;
                var dyReal = Rearth * discretization.dlat * radPerMin;
    
                var Dx = Rearth * discretization.dlon * radPerMin;
                var Dy = Rearth * discretization.dlat * radPerMin;
                var dt = Math.cos(latMax*radPerDeg) * Dy/Math.sqrt(9.81 * hmax);
    
                // let hmax = Math.max.apply(null, )
                // discretization.dt = 0.25 * Math.min(dxReal, dyReal) / Math.sqrt(9.81 * hmax);
                discretization.dt = cfl*dt;
    
            }
        }

        

        if(discretization.dt>15)
            discretization.dt = 15;

        // let hmax = Math.max.apply(null, )
        
    }

    let start = function(){

        bathymetry.texture = createTextureFromMatrix (
            bathymetry.array, bathymetry.textureId );

        if(initialSurface != undefined){
            
            initialSurface.texture = createTextureFromMatrix (
                initialSurface.array, initialSurface.textureId );

        }
        
        setTimeStep(data);

        createShaders();

        createBuffers();

        initFBOs();
        
        setPOIs();

    }

    start();

    return {
        domain,
        discretization,
        bathymetry,
        getCellBathymetry: (i,j) =>{
            return exportBuffer(wave.first.fbo, i, j, 1, 1);
        },
        getRectangleBathymetry: (iStart, jStart, iEnd, jEnd)=>{
            let rectangle = exportBuffer(wave.first.fbo, 3, iStart, jStart, iEnd-iStart+1, jEnd - jStart+1);
            rectangle = [... rectangle];
            let matrixRectangle = [];

            while (rectangle.length > 0) matrixRectangle.push(rectangle.splice(0, iEnd-iStart+1));

            return matrixRectangle;
        },
        get currentTime(){
            return discretization.stepNumber * discretization.dt;
        },
        setTimeStep,
        canvas,
        get currentGridHeights(){
            return exportBuffer(wave.first.fbo);
        },
        get currentMaximumHeights(){
            return exportBuffer(maxHeights.first.fbo);
        },
        get currentArrivalTimes(){
            return exportBuffer(maxHeights.first.fbo, 1)
        },
        pois,
        runSimulationStep,
        displayPColor: ()=>{renderDisplayProgram()},
        displayOption,
        colors
    }

}