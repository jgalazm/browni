const sphericalSWMomentumStepShaderString = `
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
              
              float hiPlusHalfj = 0.5*(uij.a + uipj.a);
              float hijPlusHalf = 0.5*(uij.a + uijp.a);
              vec2 UV = correctedUV(vUv);
              float lon = UV.x*(ymax-ymin)+ymin;
              float latj = UV.y*(ymax-ymin)+ymin;
              float coslatj = cos(degToRad(latj));
              

              float eta2ij = uij.r;
              float eta2ipj = uipj.r;
              float eta2ijp = uijp.r;


              float M2ij = 0.0;
              if(uij.a * uipj.a > gx*gx){
                  M2ij = uij.g - dt*g*hiPlusHalfj/(Rearth*coslatj*minToRad(dlon))*(eta2ipj - eta2ij);

                  // add coriolis
                  float Nc = 0.25*(uij.b + uijp.b + uipj.b + uipjp.b);
                  float R3 = 2.0 * dt * omega * sin(degToRad(latj+0.5*dlat/60.0));
                  
                  M2ij = M2ij + R3*Nc;
              }

              float N2ij = 0.0;                
              if(uij.a * uijp.a > gx*gx){
                  N2ij = uij.b - dt*g*hijPlusHalf/(Rearth*minToRad(dlat))*(eta2ijp - eta2ij);            

                  // add coriolis

                  float Mc = 0.25*(uij.g + uijp.g + uipj.g + uipjp.g);
                  float R5 = -2.0 * dt * omega * sin(degToRad(latj));
                  
                  N2ij = N2ij + R5*Mc;                            
                                      
              }

              gl_FragColor = vec4( eta2ij, M2ij, N2ij, uij.a);
          }        
      `;
export default sphericalSWMomentumStepShaderString;
