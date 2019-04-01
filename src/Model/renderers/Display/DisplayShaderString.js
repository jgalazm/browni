const DisplayShaderString = `
          precision highp float;
          
          const int nColors = 16;            
          uniform sampler2D field;
          uniform vec4 colormap[16];
          uniform float thresholds[16];
          uniform int displayedChannel;
          
          varying vec2 vUv;

          vec4 getPseudoColor(float value){
              vec4 pseudoColor;

              if(value <= thresholds[0]){
                  pseudoColor = colormap[0];
              }
              else if (value > thresholds[16-1]){
                  pseudoColor = colormap[16-1];
              }
              else{
                  for (int i=1; i<16; i++){
                      vec4 cleft = colormap[i-1];
                      vec4 cright = colormap[i];
          
                      if (value>thresholds[i-1] && value <= thresholds[i]){
                          float t = (value - thresholds[i-1])/(thresholds[i] - thresholds[i-1]);
                          pseudoColor = mix(cleft, cright, t);
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

              float h = texture2D(field, vUv).a;

              vec4 color = getPseudoColor(uij);

              color.a = color.a * step(0.0, h);

              gl_FragColor  = color;
          }    
      `;
export default DisplayShaderString;