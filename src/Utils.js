export const genericVertexShaderString = `
    precision highp  float;
    attribute vec2 inPosition;

    varying  vec2 vUv;

    void main()
    {
        vUv = inPosition.xy*0.5+0.5;
        
        gl_Position = vec4(inPosition,0, 1);
    }    
`;

export function createShaderProgram(gl, vertexShader, fragmentShader) {
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
}

export function compileShader(gl, type, source) {
  let shader = gl.createShader(type); //shader handle
  gl.shaderSource(shader, source);
  gl.compileShader(shader);

  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    console.warn(source);
    throw gl.getShaderInfoLog(shader);
  }

  return shader;
}

export function renderFrameBuffer(gl, frameBuffer) {
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuffer);
  gl.drawElements(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0);
}
