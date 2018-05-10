function Mw2Mo(Mw){
    return Math.pow(10., 1.5*Mw+9.1)
  }


function getLengthWidthSlip(Mw){

    var mu = 5.0e10
    var Mo = Mw2Mo(Mw);
    var S = 1.3E-10 *Math.pow(Mo,2/3);
    var L = Math.sqrt(2*S)*1000;
    var W = L/2.0;
    var slip = 1.66E-7 * Math.pow(Mo,1/3);
    return {L:L, W:W, slip:slip}
}

export {Mw2Mo, getLengthWidthSlip};
