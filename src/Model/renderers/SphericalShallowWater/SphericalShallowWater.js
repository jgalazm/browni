import SphericalSWMassStep from "./SphericalSWMassStep/SphericalSWMassStep";
import SphericalSWMomentumStep from "./SphericalSWMomentumStep/SphericalSWMomentumStep";

export default function SphericalShallowWater(gl) {
  const massStep = new SphericalSWMassStep(gl);

  const momentumStep = new SphericalSWMomentumStep(gl);

  const render = (doubleFBO, modelState) => {
    massStep.render(doubleFBO, modelState);

    momentumStep.render(doubleFBO, modelState);
  };

  return {
    render
  };
}
