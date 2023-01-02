
rooms.raytrace2 = function() {

lib3D();

description = `<b>Blue Moon</b>
<p><b>Change light particles</b> </p>
                  <br> <input type=range id=radius   value= 10. min="0." max="20."> light particles size
                  <br> <input type=range id=distance   value= 20. min="0." max="35."> distance of moon
                  <br>
                  <br><b>Motivation</b> : Quote from one of the books I liked.
                  <br>"That year I was 17.When I went to the wild at night, I saw that the night sky was like a pool of purple water. 
                  The stars were bright spots that did not move. The night wind is some light blue streamlines, and echos came from the clouds.
                  I was very happy at that moment, which shows that I can be a poet. 
                  According to my opinion, anyone who can feel a moment of joy in this world of endless troubles and hatreds can be regarded as a poet."

`;

code = {
'init':`





























   S.nS = 8;
   S.nL = 8;
  

  S.material = [];
  
 
   
   let materials = [
      
      [1.,1.,1.,1., 1.,1.,1.,1., 1.,1.,1.,1.,  0,0,0,0],//white
      [.6,.6,2.9,.7,     .1,.1,.1,0,  .1,5.,5.,5.,    0,0,0,0],//blue
   //   [.25,.75,0,0, 0,.8,0,0, .2,.2,2,20, 0,0,0,0],
    //  [.65,.15,.025,0, .8,.3,.05,0, 1,.6,.1,6, 0,0,0,0],
      
   ];
  
   for (let n = 0 ; n < S.nS ; n++) {
      S.material.push(materials[n % materials.length]);
   }
   
   S.velocity = [];
   S.sPos = [];
   for (let i = 0; i < S.nS; i++) {
      S.sPos.push([-.5+(i+1)*Math.random(), -1., -.5+(i+1)*Math.random()]);
      S.velocity.push([.08,.08,.08])
   }

   S.sPos2 = [];
   for (let i = 0; i < S.nS; i++) {
      S.sPos2.push([-.5+Math.random()*(i+1), .15+Math.random(), -.5+Math.random()*(i+1)]);
   }

   S.sPos3 = [];
   for (let i = 0; i < S.nS; i++) {
      S.sPos3.push([-.6+Math.random()*(i+1), .1+Math.random(), -.6+Math.random()*(i+1)]);
   }

   S.sPos4 = [];
   for (let i = 0; i < S.nS; i++) {
      S.sPos4.push([-.6+Math.random()*(i+1), .19+Math.random(), -.4+Math.random()*(i+1)]);
   }

   S.sPos5 = [];
   for (let i = 0; i < S.nS; i++) {
      S.sPos5.push([-.5, .15, -1.7]);//Position 5
   }

   S.sPos6 = [];
   for (let i = 0; i < S.nS; i++) {
      S.sPos6.push([.9, .2, -1.5]);//Position 6
   }

`,
fragment: `
S.setFragmentShader(\`
   const int nS = \` + S.nS + \`;
   const int nS2 = \` + S.nS + \`;
   const int nS3 = \` + S.nS + \`;
   const int nS4 = \` + S.nS + \`;
   const int nS5 = \` + S.nS + \`;
   const int nS6 = \` + S.nS + \`;
   const int nL = \` + S.nL + \`;
   uniform float uTime;
   uniform vec4 uS[nS];

   uniform vec4 uS2[nS2];
   uniform vec4 uS3[nS3];
   uniform vec4 uS4[nS4];
   uniform vec4 uS5[nS5];
   uniform vec4 uS6[nS6];

   uniform vec3 uLd[2];
   uniform vec3 uLc[2];
   uniform mat4 uSm[nS];
   uniform mat4 uSm2[nS2];
   uniform mat4 uSm3[nS3];
   uniform mat4 uSm4[nS4];
   uniform mat4 uSm5[nS5];
   uniform mat4 uSm6[nS6];
   uniform vec3 uDc;

   uniform mat4 uIM;
   uniform vec4 uOctahedron[8];

   varying vec3 vPos;



   float fl = 4.5;

   vec3 ld = normalize(vec3(sin(uTime),1.,1.));

   float raySphere(vec3 V, vec3 W, vec4 S) {
      V -= S.xyz;
      float b = dot(V, W);
      float d = b * b - dot(V, V) + S.w * S.w;
      return d < 0. ? -1. : -b - sqrt(d);
   }
///////////////////////////////////////////////////Begin
   // TRACE A RAY TO A HALFSPACE

   float rayHalfspace(vec3 V, vec3 W, vec4 H) {
      vec4 V1 = vec4(V, 1.);
      vec4 W0 = vec4(W, 0.);
      return -dot(H, V1) / dot(H, W0);
   }
   vec4 rayOctahedron(vec3 V, vec3 W, mat4 IM) {
      vec3 N = vec3(0.);
      float tIn = -1000., tOut = 1000.;
      for (int i = 0 ; i < 8 ; i++) {
         vec4 H = uOctahedron[i] * IM;
	 H /= sqrt(dot(H.xyz, H.xyz));
	 float t = rayHalfspace(V, W, H);
	 if (dot(W, H.xyz) < 0.) {
	    if (t > tIn)
	       N = H.xyz;
	    tIn = max(tIn, t);
	 }
	 else
	    tOut = min(tOut, t);
      }
      return vec4(N, tIn < tOut ? tIn : -1.);
   }
/////////////////////////////////////////////End 

   vec3 shadeSphere(vec3 P, vec4 S, mat4 M) {
      vec3 ambient = M[0].xyz;
      vec3 diffuse = M[0].xyz;
      vec3 specular = M[0].xyz;

      vec3 N = normalize(P - S.xyz);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {
         float t = -1.;
         for (int n = 0 ; n < nS ; n++) {
            t = max(t, raySphere(P, uLd[l], uS[n]));
         }
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular * pow(max(0., dot(R, E)), 20.));
         }
      }
      return c;
   }

   vec3 shadeSphere2(vec3 P, vec4 S, mat4 M) {
      vec3 ambient = M[0].xyz;
      vec3 diffuse = M[0].xyz;
      vec3 specular = M[0].xyz;

      vec3 N = normalize(P - S.xyz);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {
         float t = -1.;
         for (int n = 0 ; n < nS2 ; n++) {
            t = max(t, raySphere(P, uLd[l], uS2[n]));
         }
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular * pow(max(0., dot(R, E)), 20.));
         }
      }
      return c;
   }

   vec3 shadeSphere3(vec3 P, vec4 S, mat4 M) {
      vec3 ambient = M[0].xyz;
      vec3 diffuse = M[0].xyz;
      vec3 specular = M[0].xyz;

      vec3 N = normalize(P - S.xyz);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {
         float t = -1.;
         for (int n = 0 ; n < nS3 ; n++) {
            t = max(t, raySphere(P, uLd[l], uS3[n]));
         }
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular * pow(max(0., dot(R, E)), 20.));
         }
      }
      return c;
   }

   vec3 shadeSphere4(vec3 P, vec4 S, mat4 M) {
      vec3 ambient = M[0].xyz;
      vec3 diffuse = M[0].xyz;
      vec3 specular = M[0].xyz;

      vec3 N = normalize(P - S.xyz);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {
         float t = -1.;
         for (int n = 0 ; n < nS4 ; n++) {
            t = max(t, raySphere(P, uLd[l], uS4[n]));
         }
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular * pow(max(0., dot(R, E)), 20.));
         }
      }
      return c;
   }

   vec3 shadeSphere5(vec3 P, vec4 S, mat4 M) {
      vec3 ambient = M[1].xyz;
      vec3 diffuse = M[1].xyz;
      vec3 specular = M[1].xyz;

      vec3 N = normalize(P - S.xyz);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {
         float t = -1.;
         for (int n = 0 ; n < nS5 ; n++) {
            t = max(t, raySphere(P, uLd[l], uS5[n]));
         }
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular * pow(max(0., dot(R, E)), 20.));
         }
      }
      return c;
   }

   vec3 shadeSphere6(vec3 P, vec4 S, mat4 M) {
      vec3 ambient = M[1].xyz;
      vec3 diffuse = M[1].xyz;
      vec3 specular = M[1].xyz;

      vec3 N = normalize(P - S.xyz);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {
         float t = -1.;
         for (int n = 0 ; n < nS6 ; n++) {
            t = max(t, raySphere(P, uLd[l], uS6[n]));
         }
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular * pow(max(0., dot(R, E)), 20.));
         }
      }
      return c;
   }

   float turbulence(vec3 p) {
      float t = 0., f = 1.;
      for (int i = 0 ; i < 10 ; i++) {
         t += abs(noise(f * p)) / f;
         f *= 2.;
      }
      return t;
   }
   
   vec3 clouds(float y) {
      vec3 sky = .01 * vec3(.2,.0,30.9)-0.14;//color 
      float s = mix(.6,1., clamp(10.*y-2., 0.,1.));
      //sky = mix(sky, vec3(0.6,0.,0.3), min(y-1.,1.));
      sky = mix(sky, vec3(1.,0.,0.3), clamp(.1*y,0.,1.));
      sky = mix(sky, vec3(s), clamp(.5*y,0.,1.));
      return sky;

   }

   void main() {
      vec3 p = -vPos + vec3(0.02*uTime,-.19,.04*uTime);//position of cloud
      vec3 color = vec3(.17,.17,.18);
      vec3 color2 = vec3(.1,.1,.1);
      vec3 V = vec3(0.,0.,fl);
      vec3 W = normalize(vec3(vPos.xy, -fl));
      float rtMin = 10000.;
      for (int n = 0 ; n < nS ; n++) {
         float t = raySphere(V, W, uS[n]);
         if (t > 0. && t < rtMin) {
            vec3 P = V+t*W;
            color = shadeSphere(V + t * W, uS[n], uSm[n]);
            rtMin = t;

            vec3 N = normalize(P - uS[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;
         }
      }
      for (int n = 0 ; n < nS2 ; n++) {
         float t = raySphere(V, W, uS2[n]);
         if (t > 0. && t < rtMin) {
            vec3 P = V+t*W;
            color = shadeSphere2(V + t * W, uS2[n], uSm2[n]);
            rtMin = t;

            vec3 N = normalize(P - uS2[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;

         }
      }
      for (int n = 0 ; n < nS3 ; n++) {
         float t = raySphere(V, W, uS3[n]);
         if (t > 0. && t < rtMin) {
            vec3 P = V+t*W;
            color = shadeSphere2(V + t * W, uS3[n], uSm3[n]);
            rtMin = t;

            vec3 N = normalize(P - uS3[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;

         }
      }
      for (int n = 0 ; n < nS4 ; n++) {
         float t = raySphere(V, W, uS4[n]);
         if (t > 0. && t < rtMin) {
            vec3 P = V+t*W;
            color = shadeSphere2(V + t * W, uS4[n], uSm4[n]);
            rtMin = t;

            vec3 N = normalize(P - uS4[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;

         }
      }
      for (int n = 0 ; n < nS5 ; n++) {
         float t = raySphere(V, W, uS5[n]);
         if (t > 0. && t < rtMin) {
            vec3 P = V+t*W;
            color = shadeSphere2(V + t * W, uS5[n], uSm5[n]);
            rtMin = t;

            vec3 N = normalize(P - uS5[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;

         }
      }
      for (int n = 0 ; n < nS6 ; n++) {
         float t = raySphere(V, W, uS6[n]);
         if (t > 0. && t < rtMin) {
            vec3 P = V+t*W;
            color = shadeSphere2(V + t * W, uS6[n], uSm6[0]);
            rtMin = t;

            vec3 N = normalize(P - uS6[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;

         }
      }

      /////////////////////////begin
      // RAY TRACE TO THE OCTAHEDRON

      vec4 Nt = rayOctahedron(V, W, uIM);
      if (Nt.w > 0. && Nt.w < rtMin) {
         //vec3 a = mix(vec3(.1), uBgColor, .3);
         color =  vec3(max(0., dot(Nt.xyz, vec3(.5))));
      }
      ///////////////////////

      color += clouds(p.y + turbulence(p));
      gl_FragColor = vec4(sqrt(color), 1.);
   }
\`);
`,
vertex: `
S.setVertexShader(\`

   attribute vec3 aPos;
   varying   vec3 vPos;

   void main() {
      vPos = aPos;
      gl_Position = vec4(aPos, 1.);
   }

\`)
`,
render: `

// HANDY DANDY VECTOR LIBRARY

   let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
   let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   let norm = v => Math.sqrt(dot(v,v));
   
   let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
   let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];
   let matrixMultiply = function(a, b) {
      let dst = [];
      for (let n = 0 ; n < 16 ; n++)
         dst.push( a[n&3     ] * b[n&12    ] +
                   a[n&3 |  4] * b[n&12 | 1] +
                   a[n&3 |  8] * b[n&12 | 2] +
                   a[n&3 | 12] * b[n&12 | 3] );
      return dst;
   }

   
   

   let normalize = v => {
      let s = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      return [ v[0]/s, v[1]/s, v[2]/s ];
   }

   let ld0 = normalize([1,1,1]);
   let ld1 = normalize([-1,-1,-1]);




   let ldData = [];
   for (let i = 0 ; i < 3 ; i++)
      ldData.push(ld0[i]);
   for (let i = 0 ; i < 30 ; i++)
      ldData.push(ld1[i]);

   S.setUniform('3fv', 'uLd', ldData);
   S.setUniform('3fv', 'uLc', [ 1,1,1, .3,.3,.1 ]);

   /*
   S.setUniform('4fv', 'uS', [ 0,0,0,.5 ]);
   S.setUniform('4fv', 'uS2', [ 0,0,0,.5 ]);
   S.setUniform('4fv', 'uS3', [ 0,0,0,.5 ]);
   S.setUniform('4fv', 'uS4', [ 0,0,0,.5 ]);
   S.setUniform('4fv', 'uS5', [ 0,0,0,.5 ]);
  // S.setUniform('4fv', 'uS6', [ 0,0,0,.5 ]);
  */
   S.setUniform('1f', 'uTime', time);

   for (let i = 0; i < S.nS; i++) {
      S.velocity[i] = [S.velocity[i][0]+.1*(Math.random()-.5), S.velocity[i][1]+.1*(Math.random()-.8), S.velocity[i][2]+.1*(Math.random()-.5)];
      S.sPos[i][1] += .005*(i+1);
      S.sPos[i][0] += S.velocity[i][0] * .02;
      S.sPos[i][2] += S.velocity[i][2] * .02;
      if (S.sPos[i][1] > 2.) {
         S.sPos[i] = [-.7+Math.random()*(i+1)*.1, -1., -.8+Math.random()*(i+1)];
      }
   }

   for (let i = 0; i < S.nS2; i++) {
      if (S.sPos2[i][1] > 2.) {
         S.sPos2[i] = [.7*(i+1)*.1, 1., 9.];
      }
   }

   for (let i = 0; i < S.nS3; i++) {
      if (S.sPos3[i][1] > 2.) {
         S.sPos3[i] = [.6*(i+1)*.1, .7, 10.];
      }
   }

   for (let i = 0; i < S.nS4; i++) {
      if (S.sPos4[i][1] > 2.) {
         S.sPos4[i] = [.8*(i+1)*.1, .7, 10.];
      }
   }

   for (let i = 0; i < S.nS5; i++) {
      if (S.sPos5[i][1] > 2.) {
         S.sPos5[i] = [.6, 4.7, 10.];
      }
   }

   /*
   S.sVel[n][1] += .003 * Math.cos(time + (2+i) * n);
   S.sVel[n][1] += .01 * (Math.random() - .5);
   S.sPos[n][1] += .1 * S.sVel[n][i];
   */

   S.sPos6[0] = [.1*Math.sin(time), .5*Math.sin(time), -2.9];
      

   let sData = [];
   let k = 0;
   for (let n = 0 ; n < S.nS ; n++) {
      sData.push(S.sPos[k++]);
      sData.push(radius.value * Math.random()*0.0001);
   }

   let sData2 = [];
   let k2 = 0;
   for (let n = 10 ; n < 40 ; n++) {
      sData2.push(S.sPos2[k2++]);
      sData2.push(.008 * Math.sin(time));
   }

   let sData3 = [];
   let k3 = 0;
   for (let n = 10 ; n < 40 ; n++) {
      sData3.push(S.sPos3[k3++]);
      sData3.push(.008 * Math.cos(time));
   }

   let sData4 = [];
   let k4 = 0;
   for (let n = 10 ; n < 40 ; n++) {
      sData4.push(S.sPos4[k4++]);
      sData4.push(.008 * Math.sin(time)+0.0001);
   }

   let sData5 = [];
   let k5 = 0;
   for (let n = 10 ; n < 40 ; n++) {
      sData5.push(S.sPos5[k5++]);
      sData5.push(.65 *distance.value*0.01);
   }

   
   let sData6 = [];
   let k6 = 0;
   sData6.push(S.sPos6[k6++]);
   //sData6.push(.05 );
   

   S.setUniform('4fv', 'uS', sData.flat());
   S.setUniform('4fv', 'uS2', sData2.flat());
   S.setUniform('4fv', 'uS3', sData3.flat());
   S.setUniform('4fv', 'uS4', sData4.flat());
   S.setUniform('4fv', 'uS5', sData5.flat());
   S.setUniform('4fv', 'uS6', sData6.flat());
   S.setUniform('Matrix4fv', 'uSm', false, S.material.flat());
   S.setUniform('Matrix4fv', 'uSm2', false, S.material.flat());
   S.setUniform('Matrix4fv', 'uSm3', false, S.material.flat());
   S.setUniform('Matrix4fv', 'uSm4', false, S.material.flat());
   S.setUniform('Matrix4fv', 'uSm5', false, S.material[1].flat());
   S.setUniform('Matrix4fv', 'uSm6', false, S.material[1].flat());



  /////////////////////////////////////////////////////////////////
  // BEGINNINGS OF A MATRIX LIBRARY

   let matrixIdentity = () => [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];
   let matrixTranslate = (x,y,z) => [1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1];
   let matrixScale = (x,y,z) => [x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1];
   let matrixRotx = t => {
      let c = Math.cos(t), s = Math.sin(t);
      return [1,0,0,0, 0,c,s,0, 0,-s,c,0, 0,0,0,1];
   }
   let matrixRoty = t => {
      let c = Math.cos(t), s = Math.sin(t);
      return [c,0,-s,0, 0,1,0,0, s,0,c,0, 0,0,0,1];
   }
   let matrixRotz = t => {
      let c = Math.cos(t), s = Math.sin(t);
      return [c,s,0,0, -s,c,0,0, 0,0,1,0, 0,0,0,1];
   }

   let matrixInverse = function(src) {
      let dst = [], det = 0, cofactor = (c, r) => {
         let s = (i, j) => src[c+i & 3 | (r+j & 3) << 2];
         return (c+r & 1 ? -1 : 1) * ( (s(1,1) * (s(2,2) * s(3,3) - s(3,2) * s(2,3)))
                                     - (s(2,1) * (s(1,2) * s(3,3) - s(3,2) * s(1,3)))
                                     + (s(3,1) * (s(1,2) * s(2,3) - s(2,2) * s(1,3))) );
      }
      for (let n = 0 ; n < 16 ; n++) dst.push(cofactor(n >> 2, n & 3));
      for (let n = 0 ; n <  4 ; n++) det += src[n] * dst[n << 2];
      for (let n = 0 ; n < 16 ; n++) dst[n] /= det;
      return dst;
    }

    /*
   // RENDER THE POLYHEDRON

   cM = matrixMultiply(cM, matrixTranslate(Math.cos(time)/2,Math.sin(time)/2,.5));
   cM = matrixMultiply(cM, matrixRotx(time));
   cM = matrixMultiply(cM, matrixRotz(time));
   cM = matrixMultiply(cM, matrixRoty(time));
   cM = matrixMultiply(cM, matrixScale(.03,.03,.03));

   S.setUniform('Matrix4fv', 'uIM', false, matrixInverse(cM));

   
   S.setUniform('4fv', 'uOctahedron', [
      -1,-1,-1,-1, 1,1,1,-1,
      1,-1,-1,-1, -1,1,1,-1,
      -1,1,-1,-1, 1,-1,1,-1,
      -1,-1,1,-1, 1,1,-1,-1,
   ]);
   
   ////////////////////////////////////////////////////////
*/

   S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);


`,
events: `
   ;
`
};

}


