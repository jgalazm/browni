import {
  createShaderProgram,
  genericVertexShaderString,
  compileShader,
  renderFrameBuffer
} from "../../../Utils";

import MaxHeightsShaderString from "./MaxHeightsShaderString";

export default function MaxHeights(gl) {
  const maxHeightsShader = compileShader(
    gl,
    gl.FRAGMENT_SHADER,
    MaxHeightsShaderString
  );

  const genericVertexShader = compileShader(
    gl,
    gl.VERTEX_SHADER,
    genericVertexShaderString
  );

  const maxHeightsProgram = createShaderProgram(
    gl,
    genericVertexShader,
    maxHeightsShader
  );

  const render = (maxHeightsDoubleFBO, modelDoubleFBO, modelState) => {
    const { discretization } = modelState;

    gl.useProgram(maxHeightsProgram.program);
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.uniform1f(
      maxHeightsProgram.uniforms.currentTime,
      discretization.stepNumber * discretization.dt
    );
    gl.uniform1i(
      maxHeightsProgram.uniforms.maxHeights,
      maxHeightsDoubleFBO.first.textureId
    );
    gl.uniform1i(maxHeightsProgram.uniforms.wave, modelDoubleFBO.first.textureId);
    renderFrameBuffer(gl, maxHeightsDoubleFBO.second.fbo);
    maxHeightsDoubleFBO.swap();
  };

  return {
    render
  };
}
