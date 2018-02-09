var shadersCode = function (name) {
	var tsunamiVertexShaderText = `
varying vec2 vUv;
void main()
{
	vUv = uv;
	gl_Position = projectionMatrix* modelViewMatrix*vec4(position,1.0);
}
    `;

	var tsunamiInitialFragmentShaderText = `
varying vec2 vUv;
uniform sampler2D tSource;
uniform sampler2D tBati;
uniform vec2 delta;
uniform vec4 colors[7];

uniform float xmin;
uniform float xmax;
uniform float ymin;
uniform float ymax;
uniform float zmin;
uniform float zmax;
uniform float dx;
uniform float dy;

uniform float L;
uniform float W;
uniform float depth;
uniform float slip;
uniform float strike;
uniform float dip;
uniform float rake;
uniform float U3;
uniform float cn;
uniform float ce;

float I4(float db, float eta, float q, float dip, float nu, float R){
	float I = 0.0;
	if (cos(dip)>1e-24){
		I = (1.0-2.0*nu)*1.0/cos(dip)*(log(R+db)-sin(dip)*log(R+eta));
	}
	else{
		I = -(1.0-2.0*nu)*q/(R+db);
	}
	return I;

}

float I5(float xi, float eta, float q, float dip,
	float nu, float R, float db){
	float x = sqrt(xi*xi+q*q);
	float I = 0.0;
	if (cos(dip)>1e-24){
		I = (1.0-2.0*nu)*2.0/cos(dip)*
			atan(eta*(x+q*cos(dip)) + x*(R+x)*sin(dip) )
				/(xi*(R+x)*cos(dip));
	}
	else{
		I = (-1.0-2.0*nu)*xi*sin(dip)/(R+db);
	}
	return I;
}

float uz_ss(float xi, float eta, float q, float dip, float nu){
	float R = sqrt(xi*xi+eta*eta+q*q);
	float db = eta*sin(dip) - q*cos(dip);
	float u = db*q/(R*(R+eta))
			+ q*sin(dip)/(R+eta)
			+ I4(db,eta,q,dip,nu,R)*sin(dip);
	return u;
}


float uz_ds(float xi, float eta, float q, float dip, float nu){
	float R = sqrt(xi*xi+eta*eta+q*q);
	float db = eta*sin(dip) - q*cos(dip);
	float u = db*q/(R*(R+xi))
		- I5(xi,eta,q,dip,nu,R,db)*sin(dip)*cos(dip);
	if (q!=0.0){
		u = u + sin(dip)*atan(xi*eta/(q*R));
	}
	return u;
}

float uz_tf(float xi, float eta, float q, float dip, float nu){
	float R = sqrt(xi*xi+eta*eta+q*q);
	float db = eta*sin(dip) - q*cos(dip);
	float u = (eta*cos(dip)+q*sin(dip))*q/(R*(R+xi))
			+ cos(dip)*xi*q/(R*(R+eta))
			- I5(xi,eta,q,dip,nu,R,db)*sin(dip)*sin(dip);
	if (q!=0.0){
		u = u - cos(dip)*atan(xi*eta/(q*R));
	}
	return u;
}

//finally
float okada(float n, float e, float depth, float strike, float dip, float L,
			float W, float rake, float slip, float U3){
	float nu = 0.25;
	float pi = 3.14159265359;
	strike = strike*pi/180.0;
	dip = dip*pi/180.0;
	rake = rake*pi/180.0;

	float U1 = cos(rake)*slip;
	float U2 = sin(rake)*slip;

	float d = depth + sin(dip)*W/2.0;
	float ec = e ;//+ cos(strike)*cos(dip)*W/2.0;
	float nc = n ;//- sin(strike)*cos(dip)*W/2.0;
	float x = cos(strike)*nc + sin(strike)*ec +L/2.0;
	float y = sin(strike)*nc - cos(strike)*ec + cos(dip)*W/2.0;

	float p = y*cos(dip) + d*sin(dip);
	float q = y*sin(dip) - d*cos(dip);

	float F1 = uz_ss(x  ,p  ,q,dip,nu);
	float F2 = uz_ss(x  ,p-W,q,dip,nu);
	float F3 = uz_ss(x-L,p  ,q,dip,nu);
	float F4 = uz_ss(x-L,p-W,q,dip,nu);
	float F = F1-F2-F3+F4;

	float G1 = uz_ds(x  ,p  ,q,dip,nu);
	float G2 = uz_ds(x  ,p-W,q,dip,nu);
	float G3 = uz_ds(x-L,p  ,q,dip,nu);
	float G4 = uz_ds(x-L,p-W,q,dip,nu);
	float G = G1-G2-G3+G4;

	float H1 = uz_tf(x  ,p  ,q,dip,nu);
	float H2 = uz_tf(x  ,p-W,q,dip,nu);
	float H3 = uz_tf(x-L,p  ,q,dip,nu);
	float H4 = uz_tf(x-L,p-W,q,dip,nu);
	float H = H1-H2-H3+H4;


	float uz = -U1/(2.0*pi)*F - U2/(2.0*pi)*G + U3/(2.0*pi)*H;
	return uz;
}

// Stereographic projection function

vec2 stereographic_projection(float latin, float lonin,
								float lat0, float lon0){
    /*
        Gives the stereographic projection of points (lat,lon)
        onto the plane tangent to the Earth's ellipsoid in (lat0,lon0)
        Input:
            latin,lonin : coordinates of points in degrees (array)
            lat0,lon0: coordinates of tangency point in degrees (float)
    */
    float pi = 3.141592653589793 ;
    float pole = pi/2.0 - 1e-5;
    float rad_deg = pi/180.0;

    float lat = latin*rad_deg;
    float lon = lonin*rad_deg;
    float lt0 = lat0*rad_deg;
    float ln0 = lon0*rad_deg;

    lat = max(-pole,min(pole,lat));
    lt0 = max(-pole,min(pole,lt0));
    float CS = cos(lat);
    float SN = sin(lat);
    float CS0 = cos(lt0);
    float SN0 = sin(lt0);
    float xf = 0.0; //false easting
    float yf = 0.0; //false northing

    float A = 6378137.0000;            //ELLIPSOIDAL SEMI-MAJOR AXIS
    float B = 6356752.3142;            //ELLIPSOIDAL SEMI-MINOR AXIS
    float F = 0.00335281067183;      // FLATTENING, F = (A-B)/A
    float E = 0.08181919092891;        // ECCENTRICITY, E = SQRT(2.0*F-F**2)
    float F2 = 0.00669438000426;        // F2 = E**2
    float ES = 0.00673949675659;        // 2ND ECCENTRICITY, ES = E**2/(1-E**2)

    float K0 = 0.9996;                // SCALE FACTOR

    float TMP = sqrt(1.0-F2*SN*SN);
    float TMP0 = sqrt(1.0-F2*SN0*SN0);
    float RHO0 = A*(1.0-F2)/(TMP0*TMP0*TMP0);
    float NU0 = A/TMP0;
    float R = sqrt(RHO0*NU0);
    float N = sqrt(1.0+F2*CS0*CS0*CS0*CS0/(1.0-F2));

    float S1 = (1.0+SN0)/(1.0-SN0);
    float S2 = (1.0-E*SN0)/(1.0+E*SN0);
    float W1 = pow((S1*pow(S2,E)),N);
    float SN_XI0 = (W1-1.0)/(W1+1.0);
    float C = (N+SN0)*(1.0-SN_XI0)/(N-SN0)/(1.0+SN_XI0);

    float W2 = C*W1;
    float SA = (1.0+SN)/(1.0-SN);
    float SB = (1.0-E*SN)/(1.0+E*SN);
    float W = pow(C*(SA*pow(SB,E)),N);


    float XI0 = asin((W2-1.0)/(W2+1.0));
    float LM0 = ln0;

    float LM = N*(lon-LM0)+LM0;
    float XI = asin((W-1.0)/(W+1.0));

    float BETA = 1.0 + sin(XI)* sin(XI0) + cos(XI)*cos(XI0)*cos(LM-LM0);

    float Y = yf + 2.0*R*K0*(sin(XI)*cos(XI0)
                     - cos(XI)*sin(XI0)*cos(LM-LM0))/BETA;
    float X = xf + 2.0*R*K0*cos(XI)*sin(LM-LM0)/BETA;

    vec2 XY = vec2(X,Y);
    return XY;
}

void main()
{

	float n = ymin + vUv.y*(ymax-ymin);
	float e = xmin + vUv.x*(xmax-xmin);
	// n = mod(n+90.0,180.0)-90.0;
	// e = mod(e,360.0);
	
	// float cn_mod180 = mod(cn+90.0,180.0)-90.0;
	// float ce_mod360 = mod(ce,360.0);
	float cn_mod180 = cn;
	float ce_mod360 = ce;

	vec2 pos = stereographic_projection(n,e,cn_mod180,ce_mod360);

	float value = 0.0;
	if (abs(pos.x)<L/2.0*6.0 && abs(pos.y)<L/2.0*6.0){
		value  = okada(pos.g,pos.r,depth,strike,dip,L,W,rake,slip,U3);
	}
	// value = 1.0;
	float k = 40.0;
	float Lmax = max(L,W);
	value = value*smoothstep(-Lmax/2.0*k,-Lmax/2.0*2.0,pos.x);
	value = value*smoothstep(-Lmax/2.0*k,-Lmax/2.0*2.0,-pos.x);
	value = value*smoothstep(-Lmax/2.0*k,-Lmax/2.0*2.0,pos.y);
	value = value*smoothstep(-Lmax/2.0*k,-Lmax/2.0*2.0,-pos.y);


	float bati = texture2D(tBati, vUv).r;
	float bati_real = zmin + bati*(zmax-zmin);

	value = value*step(0.0,-bati_real);
	// value = 0.0;
	// if(e<=0.0){
	// 	value = 1.0;
	// }

	gl_FragColor = vec4(value,0.0,0.0,bati);

}

    `;

	var tsunamiModelFragmentShaderText = `
varying vec2 vUv;
uniform sampler2D tSource;
uniform vec2 delta;

uniform float rad_min;
uniform float rad_deg;
uniform float gx;
uniform float dx;
uniform float dy;
uniform float RR;
uniform float RS;
uniform float g;
uniform float xmin;
uniform float xmax;
uniform float ymin;
uniform float ymax;
uniform float zmin;
uniform float zmax;

vec3 getR1R11R6(float V){
	// V = vertical coordinate in [0,1];
	float ymesh = ymin + (ymax-ymin)*V;
	float angM = ymesh*rad_deg;
	float cosM = cos(angM);
	float sinM = sin(angM);
	float R1  = RR/(cosM*dx*rad_min);
	float R11 = RR/(cosM*dy*rad_min);
	float R6 = cos((ymesh+0.5*dy/60.0)*rad_deg);
	vec3 r1r11r6 = vec3(R1,R11,R6);
	return r1r11r6;
}

vec2 getR2R4(float V, float hp, float hq){
	float ymesh = ymin + (ymax-ymin)*V;
	float angM = ymesh*rad_deg;
	float cosM = cos(angM);

	float R2 = RS*max(hp,0.0)/(cosM*dx*rad_min);
	float R4 = RS*max(hq,0.0)/(cosM*dy*rad_min);

	vec2 r2r4 = vec2(R2,R4);
	return r2r4;
}


float openBoundary(vec2 vUv, float eta2_ij, vec4 u_ij, vec4 u_ijm, vec4 u_imj, float h_ij){
	float eta = eta2_ij;

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
	float k = 2.0;
	if (vUv.y <= k*delta.y){
		eta = 0.0;
		if (h_ij>gx){
			float c = sqrt(g*h_ij);
			float z = sqrt(Nij*Nij+0.25*(Mij+Mimj)*(Mij+Mimj))/c;
			if (Nij>0.0){
				z = -z;
			}
			eta = z;

		}		
	}

	if (vUv.y >= 1.0-k*delta.y){
		eta = 0.0;
		if (h_ij>gx){
			float c = sqrt(g*h_ij);
			float z = sqrt(Nijm*Nijm+0.25*(Mij+Mimj)*(Mij+Mimj))/c;
			if (Nijm<0.0){
				z = -z;
			}
			eta = z;
		}
	}

	if(xmax-xmin<360.0-0.2){
		if (vUv.x <= k*delta.x){
			eta = 0.0;
			if (h_ij>gx){
				float c = sqrt(g*h_ij);
				float z = sqrt(Mij*Mij+0.25*(Nij+Nijm)*(Nij+Nijm))/c;
				if (Mij>0.0){
					z = -z;
				}
				eta = z;

			}
			
		}
		
		if (vUv.x >=1.0-k*delta.x){
			eta = 0.0;
			if (h_ij>gx){
				float c = sqrt(g*h_ij);
				float z = sqrt(Mimj*Mimj+0.25*(Nij+Nijm)*(Nij+Nijm))/c;
				if (Mimj<0.0){
					z = -z;
				}
				eta = z;

			}
			
		}
	}
	

	if(vUv.x <= k*delta.x && vUv.y<=k*delta.y){
		eta = 0.0;
		if (h_ij>gx){
			float c = sqrt(g*h_ij);
			float z = sqrt(Mij*Mij +Nij*Nij)/c;
			if (Nij>0.0){
				z = -z;
			}
			eta = z;
		}
	}

	if(vUv.x >= 1.0-k*delta.x && vUv.y<=k*delta.y){
		eta = 0.0;
		if (h_ij>gx){
			float c = sqrt(g*h_ij);
			float z = sqrt(Mimj*Mimj+Nij*Nij)/c;
			if (Nij>0.0){
				z = -z;
			}
			eta = z;
		}
	}

	if(vUv.x <= k*delta.x && vUv.y>=1.0-k*delta.y){
		eta = 0.0;
		if (h_ij>gx){
			float c = sqrt(g*h_ij);
			float z = sqrt(Mij*Mij +Nij*Nij)/c;
			if (Nijm<0.0){
				z = -z;
			}
			eta = z;
		}
	}

	if(vUv.x >= 1.0 - k*delta.x && vUv.y>=1.0-k*delta.y){
		eta = 0.0;
		if (h_ij>gx){
			float c = sqrt(g*h_ij);
			float z = sqrt(Mimj*Mimj +Nijm*Nijm)/c;
			if (Nijm<0.0){
				z = -z;
			}
			eta = z;
		}
	}
	
	return eta;
}

float updateWaterSurface(vec2 vUv, float h_ij, vec4 u_ij, vec4 u_imj, vec4 u_ijm){
	float eta2_ij = 0.0;
	//mass equation for ij
	if (h_ij >gx){
		vec3 RRRij = getR1R11R6(vUv.y);
		float  R1ij = RRRij.x;
		float R11ij = RRRij.g;
		float R6ij = RRRij.b;

		vec3 RRRijm = getR1R11R6(vUv.y-delta.y);
		float  R1ijm = RRRijm.r;
		float R11ijm = RRRijm.g;
		float R6ijm = RRRijm.b;

		eta2_ij = u_ij.r -R1ij*(u_ij.g-u_imj.g) - R11ij*(R6ij*u_ij.b - R6ijm*u_ijm.b);
	}

	// float cn_mod180 = mod(cn+90.0,180.0)-90.0;
	// float ce_mod360 = mod(ce,360.0);
	// 
		// boundary if necessary
		eta2_ij = openBoundary(vUv, eta2_ij, u_ij, u_ijm, u_imj, h_ij);	
	
	return eta2_ij;
}

float myMod(float x){
	return x - floor(x);
}

vec2 mod1x(float x, float y){
	return vec2(myMod(x), y);
}

void main()
{
	// u = (eta, p, q, h)
	// eta: free surface
	// p: x-momentum
	// q: y-momentum
	// h: water depth, >0 if wet, <0 if dry.
	//
	// old vals
	vec4 u_ij = texture2D(tSource,vUv);

	//neighbors old vals
	vec4 u_imj = texture2D(tSource, mod1x(vUv.x - delta.x, vUv.y));
	vec4 u_ipj = texture2D(tSource, mod1x(vUv.x + delta.x, vUv.y));
	vec4 u_ijm = texture2D(tSource, mod1x(vUv.x + 0.0, vUv.y -delta.y));
	vec4 u_ijp = texture2D(tSource, mod1x(vUv.x + 0.0, vUv.y + delta.y));

	vec4 u_ipjm = texture2D(tSource, mod1x(vUv.x + delta.x,  vUv.y - delta.y));
	vec4 u_imjp = texture2D(tSource, mod1x(vUv.x + -delta.x, vUv.y + delta.y));


	//new vals
	vec4 u2_ij, u2_ipj, u2_ijp;

	//bati depth vals
	float h_ij = -(u_ij.a*(zmax-zmin)+zmin);
	float h_ipj = -(u_ipj.a*(zmax-zmin)+zmin);
	float h_ijp = -(u_ijp.a*(zmax-zmin)+zmin);
	float h_ijm = -(u_ijm.a*(zmax-zmin)+zmin);

	// update water surface in three points
	// this is needed for the momentum equation
	// (blame the leapfrog scheme)

	// update water surface for cell ij
	u2_ij.r = updateWaterSurface(vUv, h_ij, u_ij, u_imj, u_ijm);

	// update water surface for cell ij
	u2_ipj.r = updateWaterSurface(vUv+vec2(delta.x,0.0), h_ipj, u_ipj, u_ij, u_ipjm);

	// update water surface for cell ij
	u2_ijp.r = updateWaterSurface(vUv+vec2(0.0,delta.x), h_ijp, u_ijp, u_imjp, u_ij);


	float hpij = 0.5*(h_ij+h_ipj);
	float hqij = 0.5*(h_ij+h_ijp);
	vec2 RRij = getR2R4(vUv.y, h_ipj, h_ijp);
	float R2ij = RRij.r;
	float R4ij = RRij.g;

	// x -momentum
	if (vUv.y<1.0-delta.y && vUv.y>delta.y){
		if ((h_ij>gx) && (h_ipj>gx)){
			u2_ij.g = u_ij.g - R2ij*(u2_ipj.r-u2_ij.r);
		}
		else{
			u2_ij.g = 0.0;
		}
	}

	// y - momentum
	if (vUv.y<1.0-delta.y  && vUv.y>delta.y){
		if ( (h_ij>gx) && (h_ijp>gx)){
			u2_ij.b = u_ij.b - R4ij*(u2_ijp.r-u2_ij.r);
		}
		else{
			u2_ij.b = 0.0;
		}
	}

	u2_ij.a = u_ij.a;

	//sometimes .. it just blows up
	if (abs(u2_ij.r)>50.0){
		u2_ij.r = 0.0;
		u2_ij.g = 0.0;
		u2_ij.b = 0.0;
	}

	gl_FragColor = vec4(u2_ij);
}

    `;

	var tsunamiScreenFragmentShaderText = `
varying vec2 vUv;
uniform sampler2D tVisualization;
uniform vec2 delta;
const int nColors = 16;
uniform vec4 colormap[nColors];

uniform float zmin;
uniform float zmax;

vec3 getpcolor(float value){
    vec3 pseudoColor;
    //
    if(value <= colormap[0].a){
        pseudoColor = colormap[0].rgb;
    }
    else if (value > colormap[nColors-1].a){
        pseudoColor = colormap[nColors-1].rgb;
    }
    else{
        for (int i=1; i<nColors; i++){
            vec4 cleft = colormap[i-1];
            vec4 cright = colormap[i];

            if (value>cleft.a && value <=cright.a){
                float t = (value - cleft.a)/(cright.a - cleft.a);
                pseudoColor = mix(cleft.rgb, cright.rgb, t);
                break;
            }
        }
    }

	// return value*vec3(0.0,0.0,1.0);
    return pseudoColor;
}


void main()
{

    vec4 texval = texture2D(tVisualization, vUv);
    float value = texval.r;
    vec3 pseudoColor = getpcolor(value);

    float bati = texval.a;
    float bati_real = zmin + (zmax-zmin)*bati;

    float corrected_color =  pow(abs(value),0.5);
    if (bati_real>-0.0){
        corrected_color = 0.0;
    }
    else{
        float t = 1.0;
        corrected_color = t*corrected_color+1.0-t;
    }
    gl_FragColor = vec4(pseudoColor.r, pseudoColor.g,pseudoColor.b,corrected_color);
    // gl_FragColor = vec4(0.0, value,0.0,1.0);
}
    `;

	var tsunamiMaxHeightsFragmentShader = `
varying vec2 vUv;
uniform sampler2D tSource;
uniform sampler2D tMaxHeights;

void main()
{
	float bathymetry = texture2D(tSource, vUv).a;
	float currentHeight = texture2D(tSource, vUv).r;
	float previousMaxHeight = texture2D(tMaxHeights, vUv).r;

	gl_FragColor = vec4(previousMaxHeight, 0.0, 0.0, bathymetry);
	if(currentHeight>previousMaxHeight){
		gl_FragColor = vec4(currentHeight, 0.0, 0.0, bathymetry);
	}
}
    `;


	var tsunamiShadersCode = {
		vshader: tsunamiVertexShaderText,
		iFshader: tsunamiInitialFragmentShaderText,
		mFshader: tsunamiModelFragmentShaderText,
		sFshader: tsunamiScreenFragmentShaderText,
		maxHeightsFragmentShader: tsunamiMaxHeightsFragmentShader

	};


	var shaders = {
		'tsunami': tsunamiShadersCode
	}

	var getShader = function (name) {
		return shaders[name];
	}
	return getShader(name);
}
