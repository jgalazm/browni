import { getLengthWidthSlip } from "./Earthquake";
import Earthquake from './renderers/Earthquake/Earthquake';
import SphericalShallowWater from "./renderers/SphericalShallowWater/SphericalShallowWater";
import Display from "./renderers/Display/Display";
import MaxHeights from "./renderers/MaxHeights/MaxHeights";

let Model = function(data, output) {
  let gl, isWebGL2;
  let vertexShader,
    initialShader,
    asteroidShader,
    cartesianWaveShader,
    cartesianDispersiveMassShader,
    cartesianDispersiveMomentumShader,
    dispersiveMassStepShader,
    dispersiveMomentumStepShader;
  let initialProgram,
    asteroidProgram,
    cartesianWaveProgram,
    cartesianDispersiveMassProgram,
    cartesianDispersiveMomentumProgram,
    dispersiveMassStepProgram,
    dispersiveMomentumStepProgram;

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

  let compileShader = function(type, source) {
    let shader = gl.createShader(type); //shader handle
    gl.shaderSource(shader, source);
    gl.compileShader(shader);

    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      console.warn(source);
      throw gl.getShaderInfoLog(shader);
    }

    return shader;
  };

  let shaderProgram = function(vertexShader, fragmentShader) {
    let uniforms = {}; // to store uniforms handles
    let program = gl.createProgram(); //program handle

    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);

    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
      throw gl.getProgramInfoLog(program);
    }

    let uniformsCount = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);

    for (let i = 0; i < uniformsCount; i++) {
      let uniformName = gl.getActiveUniform(program, i).name;
      if (uniformName.includes("[0]"))
        uniformName = uniformName.replace("[0]", "");
      uniforms[uniformName] = gl.getUniformLocation(program, uniformName);
    }

    uniforms.vertexPositionAttribute = gl.getAttribLocation(
      program,
      "inPosition"
    );

    gl.enableVertexAttribArray(uniforms.vertexPositionAttribute);

    return { uniforms, program };
  };

  let createShaders = function() {
    vertexShader = compileShader(
      gl.VERTEX_SHADER,
      `
            precision highp  float;
            attribute vec2 inPosition;

            varying  vec2 vUv;

            void main()
            {
                vUv = inPosition.xy*0.5+0.5;
                
                gl_Position = vec4(inPosition,0, 1);
            }    
        `
    );

    initialShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
            precision highp float;
            
            varying vec2 vUv;
            
            uniform sampler2D bathymetry;
            uniform sampler2D initialSurface;

            uniform vec2 texel; 
            
            uniform float L;
            uniform float W;
            uniform float ce;
            uniform float cn;
            uniform float xmin;
            uniform float xmax;
            uniform float ymin;
            uniform float ymax;
            uniform bool isSubrectangle;

            uniform int coordinates;

            const int CARTESIAN = 0;
            const int SPHERICAL = 1;

            const float Rearth = 6378000.0;
            


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

                float eta = 0.0;
                if(isSubrectangle){
                    if (abs(pos.x)<L/2.0 && abs(pos.y)<W/2.0){
    
                        float u = (pos.x+L/2.0)/L;
                        float v = (pos.y+W/2.0)/W;
                        eta  = texture2D(initialSurface, vec2(u,v)).r;
                    }
                }
                else{    
                    eta  = texture2D(initialSurface, vUv).r;
                }

                float h = texture2D(bathymetry, vUv).r;
                h = max(0.0, h);

                gl_FragColor  = vec4(eta, 0.0, 0.0, h);
            }    
        `
    );

    asteroidShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
            precision highp float;
            varying vec2 vUv;
            uniform sampler2D bathymetry;
            
            uniform vec2 texel;

            uniform float xmin;
            uniform float xmax;
            uniform float ymin;
            uniform float ymax;

            uniform float R_i;
            uniform float v_i;
            uniform float rho_i;
            uniform float ce;
            uniform float cn;
            
            uniform int coordinates;


            const float pi = 3.14159265358979323;
            const float Rearth = 6378000.0;            
            const float g = 9.81;
            const float rho_w = 1.0;
            const float epsilon_tsunami = 0.15;
            const int CARTESIAN = 0;
            const int SPHERICAL = 1;

            vec2 simpleProjection(float latin, float lonin, float lat0, float lon0){
                float y = Rearth*(latin-lat0)*pi/180.0;
                float x = Rearth*cos(lat0*pi/180.0)*(lonin-lon0)*pi/180.0;

                return vec2(x,y);

            }

            void main(){
                float bathymetry = texture2D(bathymetry, vUv).r;
                bathymetry = max(0.0, bathymetry);


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
                    pos = simpleProjection(n,e,cn,ce);
                }


                float Q = pow(8.0*epsilon_tsunami*rho_i*v_i*v_i/(9.0*rho_w*g), 0.25);
                float Dc = Q*pow(R_i, 0.75);
                float dc = Dc*3.0;
                float Rc = dc/2.0;
                float Rd = Rc * sqrt(2.0);

                float r = length(pos);
                float height = 0.0;
                if( r <= Rd){
                    height = -Dc * (1.0-r*r/(Rc*Rc));
                }
                

                gl_FragColor = vec4(height, 0.0, 0.0, bathymetry);
            }

        `
    );
    cartesianWaveShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
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
        `
    );

    cartesianDispersiveMassShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
                precision highp float;
                
                varying vec2 vUv;

                uniform sampler2D u0;
                uniform vec2 texel;
                uniform float dt;
                uniform float dx;
                uniform float dy;

                const float eps = 1e-2;

                float massEquation(vec4 uij, vec4 uimj, vec4 uijm){
                    float eta = 0.0;
                    float hij = uij.a;
                    if(hij > eps){
                        eta = uij.r -dt/dx*(uij.g-uimj.g) -dt/dy*(uij.b-uijm.b);
                    }
                    return  eta;
                }

                void main(){
                    vec4 uij  = texture2D(u0, vUv );
                    vec4 uimj = texture2D(u0, vUv+vec2(-texel.x,0.0));
                    vec4 uijm = texture2D(u0, vUv+vec2(0.0,-texel.y));

                    vec4 u2ij = uij;

                    u2ij.r = massEquation(uij, uimj, uijm);

                    gl_FragColor = u2ij;

                }
            `
    );

    cartesianDispersiveMomentumShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
            precision highp float;
                
            varying vec2 vUv;

            uniform sampler2D u0;
            uniform vec2 texel;
            uniform float dt;
            uniform float dx;
            uniform float dy;

            const float g = 9.81;
            const float eps = 1e-2;
            
            void main(){ 
                vec2 right = vec2(texel.x, 0.0);
                vec2 front = vec2(0.0, texel.y);

                vec4 u2ij  = texture2D(u0, vUv );
                vec4 u2ipj  = texture2D(u0, vUv + right);
                vec4 u2ijp = texture2D(u0, vUv + front);

                vec4 u2ippj = texture2D(u0, vUv + 2.0*right);
                vec4 u2imj = texture2D(u0, vUv - right);
                vec4 u2ipjp = texture2D(u0, vUv + right + front);
                vec4 u2ipjm = texture2D(u0, vUv + right - front);
                vec4 u2ijpp = texture2D(u0, vUv + 2.0*front);
                vec4 u2ijm = texture2D(u0, vUv - front);
                vec4 u2imjp = texture2D(u0, vUv - right + front);

                float eta2ij = u2ij.r;
                float eta2ippj = u2ippj.r;
                float eta2ipj = u2ipj.r;
                float eta2imj = u2imj.r;
                float eta2ipjp = u2ipjp.r;
                float eta2ipjm = u2ipjm.r;
                
                float eta2ijpp = u2ijpp.r;
                float eta2ijp = u2ijp.r;
                float eta2ijm = u2ijm.r;
                float eta2imjp = u2imjp.r;

                float hij = u2ij.a;
                float hipj = u2ipj.a;
                float hijp = u2ijp.a;



                float alpha = (4.0*hij*hij + g*hij*dt*dt-dx*dx)/dx/dx;
                float gamma = alpha + 1.0;

                float M2ij = 0.0;
                if(hij>eps && hipj>eps){
                    
                    // if not dry, otherwise defaults to 0.0;
                    M2ij = u2ij.g - g*hij*dt/dx*(u2ipj.r - u2ij.r);
                    
                    
                    // dispersion correction
                    float hiPlusHalfj = 0.5*(hij + hipj);
                    M2ij = M2ij - alpha*dt/(12.0*dx)*g*hiPlusHalfj*(eta2ippj - 3.0*eta2ipj + 3.0*eta2ij - eta2imj);
                    M2ij = M2ij - gamma*dt/(12.0*dx)*g*hiPlusHalfj*(eta2ipjp -2.0*eta2ipj + eta2ipjm);
                    M2ij = M2ij + gamma*dt/(12.0*dx)*g*hiPlusHalfj*(eta2ijp - 2.0*eta2ij + eta2ijm);
                }
                
                float N2ij = 0.0;
                if(hij>eps && hijp>eps){
                    N2ij = u2ij.b - g*hij*dt/dy*(u2ijp.r - u2ij.r);
                    
                    float hijPlusHalf = 0.5*(hij + hijp);
                    // dispersion correction                    
                    N2ij = N2ij - alpha*dt/(12.0*dy)*g*hijPlusHalf*(eta2ijpp - 3.0*eta2ijp + 3.0*eta2ij - eta2ijm);
                    N2ij = N2ij - gamma*dt/(12.0*dy)*g*hijPlusHalf*(eta2ipjp - 2.0*eta2ijp + eta2imjp);
                    N2ij = N2ij + gamma*dt/(12.0*dy)*g*hijPlusHalf*(eta2ipj - 2.0*eta2ij + eta2imj);
                }
                
                gl_FragColor  = vec4(eta2ij, M2ij, N2ij, hij);
            }   
        `
    );

    dispersiveMassStepShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
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
            const float omega = 7.29e-5;


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
            
            void main(){
                // u = (eta, p, q, h)
                // eta: free surface
                // p: x-momentum
                // q: y-momentum
                // h: water depth, >0 if wet, <0 if dry.

                // read previous frame, handling periodic boundaries
                vec2 right = vec2(texel.x,0.0);
                vec2 front = vec2(0.0, texel.y);

                vec4 uij = texture2D(u0, vUv);
                vec4 uijm = texture2D(u0, vUv - front);
                vec4 uimj = texture2D(u0, vUv - right);
                if(isPeriodic == 1){
                    float uright = mod(vUv.x + right.x, 1.0);
                    float uleft = mod(vUv.x - right.x, 1.0);
                    
                    uimj = texture2D(u0, vec2(uleft, vUv.y));
                }
                
                
                
                float eta2ij = updateSurface(vUv, correctedUV(vUv), uij, uimj, uijm);


                gl_FragColor = vec4(eta2ij, uij.g, uij.b, uij.a);
            }
        `
    );

    dispersiveMomentumStepShader = compileShader(
      gl.FRAGMENT_SHADER,
      `
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
            const float omega = 7.29e-5;

            float degToRad(float degrees){
                // convert from degrees to radians
                return rad_deg*degrees;
            }

            float minToRad(float minutes){
                // convert from minutes to radians
                return rad_min*minutes;
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

            vec4 queryPeriodicTexture(vec2 uv){
                vec4 uij = texture2D(u0, uv);
                if(isPeriodic == 1){
                    float uright = mod(uv.x, 1.0);
                    uij = texture2D(u0, vec2(uright, uv.y));
                }

                return uij;
            }


            void main(){
                // u = (eta, p, q, h)
                // eta: free surface
                // p: x-momentum
                // q: y-momentum
                // h: water depth, >0 if wet, <0 if dry.
                

                // read mass step, handling periodic boundaries
                vec2 right = vec2(texel.x,0.0);
                vec2 front = vec2(0.0, texel.y);

                vec4 uij = texture2D(u0, vUv);
                vec4 uipj = queryPeriodicTexture( vUv + right );
                vec4 uipjp = queryPeriodicTexture( vUv+right+front );

                vec4 uijp = queryPeriodicTexture( vUv + front);

                vec4 u2ippj = queryPeriodicTexture( vUv + right*2.0);
                vec4 u2imj = queryPeriodicTexture( vUv - right);
                vec4 u2ipjp = queryPeriodicTexture( vUv + right + front);
                vec4 u2ipjm = queryPeriodicTexture( vUv + right - front);
                vec4 u2ijpp = queryPeriodicTexture( vUv + front*2.0);
                vec4 u2ijm = queryPeriodicTexture( vUv - front);
                vec4 u2imjp = queryPeriodicTexture( vUv - right + front);;

                float hiPlusHalfj = 0.5*(uij.a + uipj.a);
                float hijPlusHalf = 0.5*(uij.a + uijp.a);
                vec2 UV = correctedUV(vUv);
                float lon = UV.x*(ymax-ymin)+ymin;
                float latj = UV.y*(ymax-ymin)+ymin;
                float coslatj = cos(degToRad(latj));
                

                float eta2ij = uij.r;
                float eta2ipj = uipj.r;
                float eta2ijp = uijp.r;

                // float eta2ippj = u2ippj.r;
                // float eta2imj = u2imj.r;
                // float eta2ipjp = u2ipjp.r;
                // float eta2ipjm = u2ipjm.r;

                // float eta2ijpp = u2ijpp.r;
                // float eta2ijm = u2ijm.r;
                // float eta2imjp = u2imjp.r;


                // // dispersion parameters
                // float dx = Rearth*coslatj*minToRad(dlon);
                // float dy = Rearth*minToRad(dlat);
                // float alpha = (4.0*uij.a*uij.a + g * uij.a*dt*dt-dx*dx)/dx/dx;
                // float gamma = alpha + 1.0;


                float M2ij = 0.0;
                if(uij.a * uipj.a > gx*gx){
                    M2ij = uij.g - dt*g*hiPlusHalfj/(Rearth*coslatj*minToRad(dlon))*(eta2ipj - eta2ij);

                    // add coriolis
                    float Nc = 0.25*(uij.b + uijp.b + uipj.b + uipjp.b);
                    float R3 = 2.0 * dt * omega * sin(degToRad(latj+0.5*dlat/60.0));
                    
                    M2ij = M2ij + R3*Nc;

                    // dispersion correction

                    // M2ij = M2ij - alpha*dt/(12.0*dx)*g*hiPlusHalfj*(eta2ippj - 3.0*eta2ipj + 3.0*eta2ij - eta2imj);
                    // M2ij = M2ij - gamma*dt/(12.0*dx)*g*hiPlusHalfj*(eta2ipjp -2.0*eta2ipj + eta2ipjm);
                    // M2ij = M2ij + dt*(eta2ipj - 2.0*eta2ij + eta2imj);

                }

                float N2ij = 0.0;                
                if(uij.a * uijp.a > gx*gx){
                    N2ij = uij.b - dt*g*hijPlusHalf/(Rearth*minToRad(dlat))*(eta2ijp - eta2ij);            

                    // add coriolis

                    float Mc = 0.25*(uij.g + uijp.g + uipj.g + uipjp.g);
                    float R5 = -2.0 * dt * omega * sin(degToRad(latj));
                    
                    N2ij = N2ij + R5*Mc;              

                    // dispersion correction

                    // N2ij = N2ij - alpha*dt/(12.0*dy)*g*hijPlusHalf*(eta2ijpp - 3.0*eta2ijp + 3.0*eta2ij - eta2ijm);
                    // N2ij = N2ij - gamma*dt/(12.0*dy)*g*hijPlusHalf*(eta2ipjp - 2.0*eta2ijp + eta2imjp);
                    // N2ij = N2ij + dt*(eta2ipj - 2.0*eta2ij + eta2imj);                    
                                        
                }

                gl_FragColor = vec4( eta2ij, M2ij, N2ij, uij.a);
            }
        `
    );

    initialProgram = shaderProgram(vertexShader, initialShader);
    asteroidProgram = shaderProgram(vertexShader, asteroidShader);

    cartesianWaveProgram = shaderProgram(vertexShader, cartesianWaveShader);
    cartesianDispersiveMassProgram = shaderProgram(
      vertexShader,
      cartesianDispersiveMassShader
    );
    cartesianDispersiveMomentumProgram = shaderProgram(
      vertexShader,
      cartesianDispersiveMomentumShader
    );

    dispersiveMassStepProgram = shaderProgram(
      vertexShader,
      dispersiveMassStepShader
    );
    dispersiveMomentumStepProgram = shaderProgram(
      vertexShader,
      dispersiveMomentumStepShader
    );

  };

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

    gl.vertexAttribPointer(
      cartesianWaveProgram.uniforms.vertexPositionAttribute,
      2,
      gl.FLOAT,
      false,
      0,
      0
    );
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

  let renderFrameBuffer = function(frameBuffer) {
    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuffer);
    gl.drawElements(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0);
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

  /* Programs for initial condition */
  let renderInitialProgram = function() {
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );
    gl.useProgram(initialProgram.program);
    gl.uniform1i(
      initialProgram.uniforms.bathymetry,
      bathymetry.texture.textureId
    );
    gl.uniform1i(
      initialProgram.uniforms.initialSurface,
      initialSurface.texture.textureId
    );

    gl.uniform2f(
      initialProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );

    gl.uniform1f(initialProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(initialProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(initialProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(initialProgram.uniforms.ymax, domain.ymax);

    gl.uniform1f(
      initialProgram.uniforms.isSubrectangle,
      initialSurface.isSubrectangle
    );

    if (initialSurface.isSubrectangle) {
      gl.uniform1f(initialProgram.uniforms.L, initialSurface.L);
      gl.uniform1f(initialProgram.uniforms.W, initialSurface.W);
      gl.uniform1f(initialProgram.uniforms.ce, initialSurface.ce);
      gl.uniform1f(initialProgram.uniforms.cn, initialSurface.cn);
    }

    if (domain.coordinates == "cartesian") {
      gl.uniform1i(initialProgram.uniforms.coordinates, 0);
    } else if (domain.coordinates == "spherical") {
      gl.uniform1i(initialProgram.uniforms.coordinates, 1);
    }
    renderFrameBuffer(wave.first.fbo);
  };

  let renderAsteroidProgram = function() {
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );
    gl.useProgram(asteroidProgram.program);
    gl.uniform1i(
      asteroidProgram.uniforms.bathymetry,
      bathymetry.texture.textureId
    );

    gl.uniform2f(
      asteroidProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );

    gl.uniform1f(asteroidProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(asteroidProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(asteroidProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(asteroidProgram.uniforms.ymax, domain.ymax);

    gl.uniform1f(asteroidProgram.uniforms.R_i, asteroid.R_i);
    gl.uniform1f(asteroidProgram.uniforms.v_i, asteroid.v_i);
    gl.uniform1f(asteroidProgram.uniforms.rho_i, asteroid.rho_i);
    gl.uniform1f(asteroidProgram.uniforms.ce, asteroid.ce);
    gl.uniform1f(asteroidProgram.uniforms.cn, asteroid.cn);

    if (domain.coordinates == "cartesian") {
      gl.uniform1i(asteroidProgram.uniforms.coordinates, 0);
    } else if (domain.coordinates == "spherical") {
      gl.uniform1i(asteroidProgram.uniforms.coordinates, 1);
    }

    renderFrameBuffer(wave.first.fbo);
  };

  /* Cartesian equations */
  let renderCartesianProgram = function() {
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );
    gl.useProgram(cartesianWaveProgram.program);
    gl.uniform2f(
      cartesianWaveProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1i(cartesianWaveProgram.uniforms.u0, wave.first.textureId);

    gl.uniform1f(cartesianWaveProgram.uniforms.dx, discretization.dx);
    gl.uniform1f(cartesianWaveProgram.uniforms.dy, discretization.dy);
    gl.uniform1f(cartesianWaveProgram.uniforms.dt, discretization.dt);

    renderFrameBuffer(wave.second.fbo);
    wave.swap();
  };

  /* Cartesian dispersive equations */
  let renderCartesianDispersiveMassProgram = function() {
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.useProgram(cartesianDispersiveMassProgram.program);

    gl.uniform1i(
      cartesianDispersiveMassProgram.uniforms.u0,
      wave.first.textureId
    );
    gl.uniform2f(
      cartesianDispersiveMassProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1f(cartesianDispersiveMassProgram.uniforms.dt, discretization.dt);
    gl.uniform1f(cartesianDispersiveMassProgram.uniforms.dx, discretization.dx);
    gl.uniform1f(cartesianDispersiveMassProgram.uniforms.dy, discretization.dy);

    renderFrameBuffer(wave.second.fbo);
    wave.swap();
  };

  let renderCartesianDispersiveMomentumProgram = function() {
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.useProgram(cartesianDispersiveMomentumProgram.program);

    gl.uniform1i(
      cartesianDispersiveMomentumProgram.uniforms.u0,
      wave.first.textureId
    );
    gl.uniform2f(
      cartesianDispersiveMomentumProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1f(
      cartesianDispersiveMomentumProgram.uniforms.dt,
      discretization.dt
    );
    gl.uniform1f(
      cartesianDispersiveMomentumProgram.uniforms.dx,
      discretization.dx
    );
    gl.uniform1f(
      cartesianDispersiveMomentumProgram.uniforms.dy,
      discretization.dy
    );

    renderFrameBuffer(wave.second.fbo);
    wave.swap();
  };

  let renderCartesianDispersiveProgram = function() {
    renderCartesianDispersiveMassProgram();
    renderCartesianDispersiveMomentumProgram();
  };

  /* Spherical dispersive equations */
  let renderDispersiveMassStepProgram = function() {
    gl.useProgram(dispersiveMassStepProgram.program);
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.uniform1f(dispersiveMassStepProgram.uniforms.dlon, discretization.dlon);
    gl.uniform1f(dispersiveMassStepProgram.uniforms.dlat, discretization.dlat);
    gl.uniform1f(dispersiveMassStepProgram.uniforms.dt, discretization.dt);
    gl.uniform1f(dispersiveMassStepProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(dispersiveMassStepProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(dispersiveMassStepProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(dispersiveMassStepProgram.uniforms.ymax, domain.ymax);

    gl.uniform2f(
      dispersiveMassStepProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1i(dispersiveMassStepProgram.uniforms.u0, wave.first.textureId);
    gl.uniform1i(
      dispersiveMassStepProgram.uniforms.isPeriodic,
      domain.isPeriodic
    );
    renderFrameBuffer(wave.second.fbo);
    wave.swap();
  };

  let renderDispersiveMomentumStepProgram = function() {
    gl.useProgram(dispersiveMomentumStepProgram.program);
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.uniform1f(
      dispersiveMomentumStepProgram.uniforms.dlon,
      discretization.dlon
    );
    gl.uniform1f(
      dispersiveMomentumStepProgram.uniforms.dlat,
      discretization.dlat
    );
    gl.uniform1f(dispersiveMomentumStepProgram.uniforms.dt, discretization.dt);
    gl.uniform1f(dispersiveMomentumStepProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(dispersiveMomentumStepProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(dispersiveMomentumStepProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(dispersiveMomentumStepProgram.uniforms.ymax, domain.ymax);

    gl.uniform2f(
      dispersiveMomentumStepProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1i(
      dispersiveMomentumStepProgram.uniforms.u0,
      wave.first.textureId
    );
    gl.uniform1i(
      dispersiveMomentumStepProgram.uniforms.isPeriodic,
      domain.isPeriodic
    );
    renderFrameBuffer(wave.second.fbo);
    wave.swap();
  };

  let renderDispersiveWaveEquation = function() {
    renderDispersiveMassStepProgram();
    renderDispersiveMomentumStepProgram();
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

      // pois[poi].shallowCorrectionFactor = 1;
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

    maxHeightsStep.render( maxHeights, wave, modelState);

    storePOISValues();
  };

  let setTimeStep = function(options) {
    let hmax = Math.max.apply(
      null,
      bathymetry.array.map(row => {
        return Math.max.apply(null, row);
      })
    );

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
    bathymetry.texture = createTextureFromMatrix(
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

    createShaders();

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
