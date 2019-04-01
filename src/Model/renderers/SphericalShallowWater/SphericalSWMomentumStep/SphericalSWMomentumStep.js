import {
  createShaderProgram,
  genericVertexShaderString,
  compileShader,
  renderFrameBuffer
} from "../../../../Utils";

import sphericalSWMomentumStepShaderString from "./SphericalSWMomentumStepShaderString";

export default function SphericalSWMomentumStep(gl) {

  const genericVertexShader = compileShader(
    gl,
    gl.VERTEX_SHADER,
    genericVertexShaderString
  );

  const sphericalSWMomentumStepShader = compileShader(
    gl,
    gl.FRAGMENT_SHADER,
    sphericalSWMomentumStepShaderString
  );

  const sphericalSWMomentumStepProgram = createShaderProgram(
    gl,
    genericVertexShader,
    sphericalSWMomentumStepShader
  );

  const render = (doubleFBO, modelState) => {
    const {discretization, domain} = modelState;
    gl.useProgram(sphericalSWMomentumStepProgram.program);
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.uniform1f(
      sphericalSWMomentumStepProgram.uniforms.dlon,
      discretization.dlon
    );
    gl.uniform1f(
      sphericalSWMomentumStepProgram.uniforms.dlat,
      discretization.dlat
    );
    gl.uniform1f(sphericalSWMomentumStepProgram.uniforms.dt, discretization.dt);
    gl.uniform1f(sphericalSWMomentumStepProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(sphericalSWMomentumStepProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(sphericalSWMomentumStepProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(sphericalSWMomentumStepProgram.uniforms.ymax, domain.ymax);

    gl.uniform2f(
      sphericalSWMomentumStepProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1i(
      sphericalSWMomentumStepProgram.uniforms.u0,
      doubleFBO.first.textureId
    );
    gl.uniform1i(
      sphericalSWMomentumStepProgram.uniforms.isPeriodic,
      domain.isPeriodic
    );
    renderFrameBuffer(gl, doubleFBO.second.fbo);
    doubleFBO.swap();
  };

  return {
    render
  };
}
