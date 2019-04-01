const MaxHeightsShaderString = `
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

              float depth = texture2D(wave, vUv).a;

              if( maxHeight < newHeight){
                  maxHeight = newHeight;
              }
              
              if(arrivalTime == 0.0 && newHeight > 1e-5 && depth > 100.0){
                  arrivalTime = currentTime;
              }
              gl_FragColor = vec4( maxHeight, arrivalTime, 0.0, 1.0);

          }
      `;
export default MaxHeightsShaderString;