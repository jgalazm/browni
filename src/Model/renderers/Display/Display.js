import {
  createShaderProgram,
  genericVertexShaderString,
  compileShader,
  renderFrameBuffer
} from "../../../Utils";

import DisplayShaderString from "./DisplayShaderString";

export default function Display(gl) {
  const displayShader = compileShader(
    gl,
    gl.FRAGMENT_SHADER,
    DisplayShaderString
  );

  const genericVertexShader = compileShader(
    gl,
    gl.VERTEX_SHADER,
    genericVertexShaderString
  );

  const displayProgram = createShaderProgram(
    gl,
    genericVertexShader,
    displayShader
  );

  const render = (doubleFBO, displayOption, colormap) => {

    gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
    gl.useProgram(displayProgram.program);
    gl.uniform4fv(
      displayProgram.uniforms.colormap,
      new Float32Array(colormap.rgba)
    );
    gl.uniform1fv(
      displayProgram.uniforms.thresholds,
      new Float32Array(colormap.thresholds)
    );

    let displayedChannel = 0;
    if (displayOption === "heights") {
      gl.uniform1i(displayProgram.uniforms.field, doubleFBO.first.textureId); //TDDO: fix texid
    } else if (
      displayOption === "max heights" ||
      displayOption === "arrival times"
    ) {
      gl.uniform1i(displayProgram.uniforms.field, maxHeights.first.textureId);
      if (displayOption === "arrival times") displayedChannel = 1;
    }
    gl.uniform1i(displayProgram.uniforms.displayedChannel, displayedChannel);

    renderFrameBuffer(gl, null);
  };

  return {
    render
  };
}