const okadaShaderString = `
            precision highp float;
            varying vec2 vUv;
            uniform sampler2D bathymetry;
            
            uniform vec2 texel; 
            
            uniform float xmin;
            uniform float xmax;
            uniform float ymin;
            uniform float ymax;
            uniform float bathymetry_xmin;
            uniform float bathymetry_ymin;
            uniform float bathymetry_xmax;
            uniform float bathymetry_ymax;

            uniform sampler2D previousTexture;
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
            uniform int reference;
            uniform int coordinates;
            
            const float MY_PI = 3.14159265358979323;
            const float Rearth = 6378000.0;
            // const float Rearth = 6384.e+3;
            const int CENTER = 0;
            const int BEGINTOP = 1;
            const int MIDTOP = 2;
            const int BEGINBOTTOM = 3;
            const int MIDBOTTOM = 4;

            const int CARTESIAN = 0;
            const int SPHERICAL = 1;

            float I4(float db, float eta, float q, float dip, float nu, float R){
                float I = 0.0;
                if (abs(cos(dip))>1e-24){
                    I = (1.0-2.0*nu)/cos(dip)*(log(R+db)-sin(dip)*log(R+eta));
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
                if (abs(cos(dip))>1e-24){
                    if(abs(xi)<1e-10){
                        I = 0.0;
                    }
                    else{
                        I = (1.0-2.0*nu)*2.0/cos(dip)*
                            atan( 
                                (  eta*(x+q*cos(dip)) + x*(R+x)*sin(dip)  )
                                /(  xi*(R+x)*cos(dip)   )
                            );                        
                    }
                
                }
                else{
                    I = -(1.0-2.0*nu)*xi*sin(dip)/(R+db);
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
                float R = sqrt(xi*xi+ eta*eta+ q*q);
                float db = eta*sin(dip) - q*cos(dip);
                float u = db*q/(R*(R+xi))
                    - I5(xi,eta,q,dip,nu,R,db)*sin(dip)*cos(dip);
                if (q*R!=0.0){
                    u = u + sin(dip)*atan(xi*eta/(q*R));
                }
                else{
                    if(xi*eta>0.0){
                        u = u + MY_PI*sin(dip);
                    }
                    else if(xi*eta<0.0){
                        u = u - MY_PI*sin(dip);
                    }
                }
                return u;
            }

            float uz_tf(float xi, float eta, float q, float dip, float nu){
                float R = sqrt(xi*xi+eta*eta+q*q);
                float db = eta*sin(dip) - q*cos(dip);
                float u = (eta*cos(dip)+q*sin(dip))*q/(R*(R+xi))
                        + cos(dip)*xi*q/(R*(R+eta))
                        - I5(xi,eta,q,dip,nu,R,db)*sin(dip)*sin(dip);
                if (q*R!=0.0){
                    u = u - cos(dip)*atan(xi*eta/(q*R));
                }
                else{
                    if(xi*eta>0.0){
                        u = u - MY_PI*cos(dip);
                    }
                    else if(xi*eta<0.0){
                        u = u + MY_PI*cos(dip);
                    }
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
                
                // const int CENTER = 0;
                // const int BEGINTOP = 1;
                // const int MIDTOP = 2;
                // const int BEGINBOTTOM = 3;
                // const int MIDBOTTOM = 4;
        
                // rotate (n,e) coordinates centered at the reference point (after stereo proj)
                float x = cos(strike)*n + sin(strike)*e;
                float y = sin(strike)*n - cos(strike)*e;
                float d = depth;

                // translate according to reference point location
                // begin bottom is the default case for this formula
                if(reference==CENTER){
                    x = x + 0.5*L;
                    y = y + 0.5*W*cos(dip);
                    d = d + 0.5*W*sin(dip);         
                }
                else if(reference == BEGINTOP){
                    y = y + W*cos(dip);
                    d = d + W*sin(dip);
                }
                else if(reference == MIDTOP){
                    x = x + 0.5*L;
                    y = y + W*cos(dip);
                    d = d + W*sin(dip);
                }
                else if(reference == MIDBOTTOM){
                    x = x + 0.5*L;
                }


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


                float uz = -U1/(2.0*pi)*F - U2/(2.0*pi)*G + U3/(2.0*pi)*H ;

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

            vec2 simpleProjection(float latin, float lonin, float lat0, float lon0){
                float pi = 3.141592653589793 ;
                float y = Rearth*(latin-lat0)*pi/180.0;
                float x = Rearth*cos(lat0*pi/180.0)*(lonin-lon0)*pi/180.0;

                return vec2(x,y);

            }

            void main()
            {

                // normalize vUv so it is defined 1-1 in [0,1] 
                // this way the first pixel corresponds to xmin and the last to xmax, exactly

                float nx = 1.0/texel.x;
                float ny = 1.0/texel.y;
                float V = (vUv.y-0.5*texel.y)/((ny-1.0)*texel.y);
                float U = (vUv.x-0.5*texel.x)/((nx-1.0)*texel.x);
                float n = ymin + V*(ymax-ymin);
                float e = xmin + U*(xmax-xmin);

                // center on reference point and make projection if necessary
                vec2 pos;
                if(coordinates==CARTESIAN){
                    pos = vec2(e-ce,n-cn);
                }
                else if(coordinates==SPHERICAL){
                    // vec2 pos = stereographic_projection(n,e,cn,ce);
                    pos = simpleProjection(n,e,cn,ce);
                }
                
                // only calculate in a square around the rectangle
                float value = 0.0;
                if (abs(pos.x)<L/2.0*8.0 && abs(pos.y)<L/2.0*8.0){

                    value  = okada(pos.g, pos.r, depth, strike, dip, L, W, rake, slip, U3);

                }

                float u_bathymetry = (e - bathymetry_xmin)/(bathymetry_xmax - bathymetry_xmin);
                float v_bathymetry = (n - bathymetry_ymin)/(bathymetry_ymax - bathymetry_ymin);

                float bathymetry = texture2D(bathymetry, vec2(u_bathymetry, v_bathymetry)).r;
                value = value*step(0.0,bathymetry);
                bathymetry = max(0.0, bathymetry);

                float u0 = texture2D(previousTexture, vUv).r;

                value = value + u0;


                gl_FragColor = vec4(value,0.0,0.0, bathymetry+value);

            }
        
        `;

export default okadaShaderString;
