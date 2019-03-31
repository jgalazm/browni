import {
  createShaderProgram,
  genericVertexShaderString,
  compileShader,
  renderFrameBuffer
} from "../../../../Utils";

import sphericalSWMassStepShaderString from "./SphericalSWMassStepShaderString";

export default function SphericalSWMassStep(gl) {

  const genericVertexShader = compileShader(
    gl,
    gl.VERTEX_SHADER,
    genericVertexShaderString
  );

  const sphericalSWMassStepShader = compileShader(
    gl,
    gl.FRAGMENT_SHADER,
    sphericalSWMassStepShaderString
  );

  const sphericalSWMassStepProgram = createShaderProgram(
    gl,
    genericVertexShader,
    sphericalSWMassStepShader
  );

  const render = (doubleFBO, modelState) => {
    const {discretization, domain} = modelState;

    gl.useProgram(sphericalSWMassStepProgram.program);
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );

    gl.uniform1f(sphericalSWMassStepProgram.uniforms.dlon, discretization.dlon);
    gl.uniform1f(sphericalSWMassStepProgram.uniforms.dlat, discretization.dlat);
    gl.uniform1f(sphericalSWMassStepProgram.uniforms.dt, discretization.dt);
    gl.uniform1f(sphericalSWMassStepProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(sphericalSWMassStepProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(sphericalSWMassStepProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(sphericalSWMassStepProgram.uniforms.ymax, domain.ymax);

    gl.uniform2f(
      sphericalSWMassStepProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );
    gl.uniform1i(sphericalSWMassStepProgram.uniforms.u0, doubleFBO.first.textureId);
    gl.uniform1i(
      sphericalSWMassStepProgram.uniforms.isPeriodic,
      domain.isPeriodic
    );
    renderFrameBuffer(gl, doubleFBO.second.fbo);
    doubleFBO.swap();
  };

  return {
    render
  };
}
