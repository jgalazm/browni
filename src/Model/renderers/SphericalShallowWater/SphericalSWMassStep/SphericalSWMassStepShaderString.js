const sphericalSWMassStepShaderString = `
precision highp float;

varying vec2 vUv;
uniform sampler2D u0; 
uniform vec2 texel; 

uniform float dlon;
uniform float dlat;
uniform float dt;
uniform float xmin; 
uniform float xmax;
uniform float ymin;
uniform float ymax;
uniform int isPeriodic;

const float rad_min = 0.000290888208665721; 
const float rad_deg = 0.01745329252;
const float cori_w = 7.2722e-5;
const float gx = 1e-5;
const float g = 9.81;
const float Rearth = 6378000.0;
const float omega = 7.29e-5;


float degToRad(float degrees){
    // convert from degrees to radians
    return rad_deg*degrees;
}

float minToRad(float minutes){
    // convert from minutes to radians
    return rad_min*minutes;
}        

float openBoundary(vec2 vUv, vec4 u_ij, vec4 u_ijm, vec4 u_imj, float h_ij){
    float eta = 0.0;
    if(h_ij<gx){
        return eta;
    }
    float c = sqrt(g*h_ij);

    float etaij = u_ij.r;
    float Mij = u_ij.g;
    float Nij = u_ij.b;
    
    float etaimj = u_imj.r;
    float Mimj = u_imj.g;
    float Nimj = u_imj.b;

    float etaijm = u_ijm.r;
    float Mijm = u_ijm.g;
    float Nijm = u_ijm.b;
    
    //j=0
    float k = 1.0;
    if (vUv.y <= k*texel.y){
        eta = sqrt(Nij*Nij+0.25*(Mij+Mimj)*(Mij+Mimj))/c;

        if (Nij>0.0){
            eta = -eta;
        }
    }

    if (vUv.y >= 1.0-k*texel.y){
        // eta = sqrt(Nijm*Nijm+0.25*(Mij+Mijm)*(Mij+Mijm))/c;
        eta = sqrt(Nijm*Nijm + Mijm*Mijm)/c;
        if (Nijm<0.0){
            eta = -eta;
        }
    }

    if(isPeriodic==0){
        if (vUv.x <= k*texel.x){
            eta = sqrt(Mij*Mij+0.25*(Nij+Nijm)*(Nij+Nijm))/c;
            if (Mij>0.0){
                eta = -eta;
            }
        }
        
        if (vUv.x >=1.0-k*texel.x){
            eta = sqrt(Mimj*Mimj+0.25*(Nij+Nijm)*(Nij+Nijm))/c;
            if (Mimj<0.0){
                eta = -eta;
            }
        }
    }
    

    if(vUv.x <= k*texel.x && vUv.y<=k*texel.y){
        eta = sqrt(Mij*Mij +Nij*Nij)/c;
        if (Nij>0.0){
            eta = -eta;
        }
    }

    if(vUv.x >= 1.0-k*texel.x && vUv.y<=k*texel.y){
        eta = sqrt(Mimj*Mimj+Nij*Nij)/c;
        if (Nij>0.0){
            eta = -eta;
        }
    }

    if(vUv.x <= k*texel.x && vUv.y>=1.0-k*texel.y){
        eta = sqrt(Mij*Mij +Nij*Nij)/c;
        if (Nijm<0.0){
            eta = -eta;
        }
    }

    if(vUv.x >= 1.0 - k*texel.x && vUv.y>=1.0-k*texel.y){
        eta = sqrt(Mimj*Mimj +Nijm*Nijm)/c;
        if (Nijm<0.0){
            eta = -eta;
        }
    }

    if(abs(eta)>1.0){
      eta = 0.0;
    }
    
    return eta;
}

float updateInnerCellSurface(vec2 vUv, vec2 UV, vec4 uij, vec4 uimj, vec4 uijm){
    /* apply mass conservation equation to update water surface */

    // calculate space varying constants
    float V  = UV.y;
    float dlonrad = minToRad(dlon);
    float dlatrad = minToRad(dlat);
    float latj = V*(ymax-ymin)+ymin;
    float latjp = latj + 0.5*dlat/60.0;
    float latjm = latj - 0.5*dlat/60.0;
    float coslatj = cos(degToRad(latj));
    float coslatjp = cos(degToRad(latjp));
    float coslatjm = cos(degToRad(latjm));

    float eta2ij = 0.0;

    if(uij.a>gx){
        eta2ij = uij.r - dt/(Rearth*coslatj*dlonrad)*(uij.g - uimj.g + uij.b*coslatjp - uijm.b*coslatjm); 
    }

    return eta2ij;

}

float updateSurface(vec2 vUv, vec2 UV, vec4 uij, vec4 uimj, vec4 uijm){

    bool isBoundary = vUv.x < texel.x || vUv.x > 1.0 - texel.x || vUv.y < texel.y || vUv.y > 1.0 - texel.y;

   
    if (isBoundary && isPeriodic == 0){
        
        return openBoundary(UV, uij, uijm, uimj, uij.a);


    }
    else{

        return updateInnerCellSurface(vUv, UV, uij, uimj, uijm);

    }               

    return uij.r;

}
vec4 queryPeriodicTexture( vec2 uv){
    vec4 uij = texture2D(u0, uv);
    if(isPeriodic == 1){
        float uright = mod(uv.x, 1.0);
        uij = texture2D(u0, vec2(uright, uv.y));
    }

    return uij;
}
vec2 correctedUV(vec2 uv){
    // corrected uv coordinates 
    // so min(u) maps to U=0, and max(u) maps to U = 1
    // same with v
    // since normally, these are displace by half texel

    float nx = 1.0/texel.x;
    float ny = 1.0/texel.y;
    float V = (uv.y-0.5*texel.y)/((ny-1.0)*texel.y);
    float U = (uv.x-0.5*texel.x)/((nx-1.0)*texel.x);

    return vec2(U,V);
}

void main(){
    // u = (eta, p, q, h)
    // eta: free surface
    // p: x-momentum
    // q: y-momentum
    // h: water depth, >0 if wet, <0 if dry.

    // read previous frame, handling periodic boundaries
    vec2 right = vec2(texel.x,0.0);
    vec2 front = vec2(0.0, texel.y);

    vec4 uij = texture2D(u0, vUv);
    vec4 uijm = texture2D(u0, vUv - front);
    vec4 uimj = queryPeriodicTexture(vUv - right);

    float eta2ij = updateSurface(vUv, correctedUV(vUv), uij, uimj, uijm);


    gl_FragColor = vec4(eta2ij, uij.g, uij.b, uij.a);
    // gl_FragColor = uimj;
}

`;

export default sphericalSWMassStepShaderString;