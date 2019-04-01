import Okada from "./Okada";

export default function Earthquake(gl) {
  const okadaRenderer = new Okada(gl);

  const render = (doubleFBO, modelState, earthquake) => {
    earthquake.forEach(finiteFault => {
      okadaRenderer.render(doubleFBO, modelState, finiteFault);
      doubleFBO.swap();
    });
  };

  return {
    render
  };
}
