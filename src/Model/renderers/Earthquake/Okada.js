import {
  createShaderProgram,
  genericVertexShaderString,
  compileShader,
  renderFrameBuffer
} from "../../../Utils";

import okadaShaderString from './OkadaShaderString';

export default function Okada(gl) {
  
  const okadaShader = compileShader(gl, gl.FRAGMENT_SHADER, okadaShaderString);

  const vertexShader = compileShader(
    gl,
    gl.VERTEX_SHADER,
    genericVertexShaderString
  );

  const okadaProgram = createShaderProgram(gl, vertexShader, okadaShader);

  const render = (doubleFBO, modelState, finiteFault) => {
    const { discretization, bathymetry, domain } = modelState;
    gl.viewport(
      0,
      0,
      discretization.numberOfCells[0],
      discretization.numberOfCells[1]
    );
    gl.useProgram(okadaProgram.program);
    gl.uniform1i(
      okadaProgram.uniforms.bathymetry,
      bathymetry.texture.textureId
    );

    gl.uniform2f(
      okadaProgram.uniforms.texel,
      1 / discretization.numberOfCells[0],
      1 / discretization.numberOfCells[1]
    );

    gl.uniform1f(okadaProgram.uniforms.xmin, domain.xmin);
    gl.uniform1f(okadaProgram.uniforms.xmax, domain.xmax);
    gl.uniform1f(okadaProgram.uniforms.ymin, domain.ymin);
    gl.uniform1f(okadaProgram.uniforms.ymax, domain.ymax);
    gl.uniform1f(okadaProgram.uniforms.bathymetry_xmin, bathymetry.extent.xmin);
    gl.uniform1f(okadaProgram.uniforms.bathymetry_ymin, bathymetry.extent.ymin);
    gl.uniform1f(okadaProgram.uniforms.bathymetry_xmax, bathymetry.extent.xmax);
    gl.uniform1f(okadaProgram.uniforms.bathymetry_ymax, bathymetry.extent.ymax);

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
    gl.uniform1i(
      okadaProgram.uniforms.previousTexture,
      doubleFBO.first.textureId
    );

    finiteFault.reference = finiteFault.reference.replace(
      String.fromCharCode(13),
      ""
    );

    if (finiteFault.reference == "center") {
      gl.uniform1i(okadaProgram.uniforms.reference, 0);
    } else if (finiteFault.reference == "begin top") {
      gl.uniform1i(okadaProgram.uniforms.reference, 1);
    } else if (finiteFault.reference == "mid top") {
      gl.uniform1i(okadaProgram.uniforms.reference, 2);
    } else if (finiteFault.reference == "begin bottom") {
      gl.uniform1i(okadaProgram.uniforms.reference, 3);
    } else if (finiteFault.reference == "mid bottom") {
      gl.uniform1i(okadaProgram.uniforms.reference, 4);
    }

    if (domain.coordinates == "cartesian") {
      gl.uniform1i(okadaProgram.uniforms.coordinates, 0);
    } else if (domain.coordinates == "spherical") {
      gl.uniform1i(okadaProgram.uniforms.coordinates, 1);
    }

    renderFrameBuffer(gl, doubleFBO.second.fbo);
  };

  return {
    render
  };
}
