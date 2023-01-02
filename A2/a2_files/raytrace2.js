
rooms.raytrace2 = function() {

lib3D();

description = `Raytracing to spheres<br>in a fragment shader


<br> <b>Background color</b>
    <br> <input type=range id=red   value= 0> red
    <br> <input type=range id=green value=0> green
    <br> <input type=range id=blue  value=0> blue</br>


<br><b>change light color</b>
<br><input type=range id=lred  value=100 > light red
    <br> <input type=range id=lgreen > light green
    <br> <input type=range id=lblue  > light blue

`;

code = {
'init':`
   S.nS = 45;//sphere number
   S.nL = 3; //light number

 let materials = [
   //   [.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,3, 0,0,0,0], // COPPER
   //   [.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0], // GOLD
   //   [.25,0,0,0,      .5,0,0,0,    2,2,2,20,   0,0,0,0], // PLASTIC
   //   [.05,.05,.05,0,  .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // LEAD
   //   [.1,.1,.1,0,     .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // SILVER
      [.1,.3,.6,0,     .1,.4,.2,0,  1,1,1,5,    0,0,0,0], // blue 
[.2,.3,.1,0,     .1,.4,.2,0,  5,5,5,30,    0,0,0,0], // green
[.1,.3,.6,0,     1,.6,.1,6,  1,1,1,5,    0,0,0,0], // white 
   ];

//make sure all spheres have material
   S.material = [];
   for (let n = 0 ; n < S.nS ; n++)
      S.material.push(materials[n % materials.length]);

/*
   S.sPos = [];
   for (let n = 0 ; n < S.nS ; n++)
      for (let i = 0 ; i < 10 ; i++)
         S.sPos.push(Math.random() - .5);
*/
   S.sPos = [];
   S.sVel = [];
   for (let n = 0 ; n < S.nS ; n++) {
      S.sPos.push([ Math.random() - .5,
                    Math.random() - .5,
                    Math.random() - .5 ]);
      S.sVel.push([0,0,0]);
   }


`,
fragment: `
S.setFragmentShader(\`
   const int nS = \` + S.nS + \`;
   const int nL = \` + S.nL + \`;



   uniform float uTime;
   uniform vec3 uBgColor;
   uniform vec4 uS[nS];
   uniform mat4 uSm[nS];
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];

   varying vec3 vPos;

   float fl = 3.;

   vec3 ld = normalize(vec3(sin(10.*uTime),1.,1.));

   float raySphere(vec3 V, vec3 W, vec4 S) {
      V -= S.xyz;
	V+=.01*W;
      float b = dot(V, W);
      float d = b * b - dot(V, V) + S.w * S.w;
      return d < 0. ? -1. : -b - sqrt(d);
   }

   vec3 shadeSphere(vec3 P, vec4 S, mat4 M) {

vec3 ambient =M[0].rgb;
vec3 diffuse =M[1].rgb;
vec3 specular =M[2].rgb;
float p =M[2].a;

      vec3 N = normalize(P - S.xyz);
      //vec3 c = vec3(.1,.1,.2);

 vec3 c = mix(ambient, uBgColor, .3);
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {

         float t = -1.;
	 for (int n = 0 ; n < nS ; n++)
	    t = max(t, raySphere(P, uLd[l], uS[n]));

         float weight = step(t, 0.);
         if (t < 0.) {
//highlight
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += weight * uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
	                   + specular*pow(max(0., dot(R, E)), 20.));
         }
      }
      //c *= 1. + .5 * noise(3.*N);
      return c;
   }

   void main() {
      vec3 color = uBgColor;
      vec3 V = vec3(0.,0.,fl);
      vec3 W = normalize(vec3(vPos.xy, -fl));

      float tMin = 10000.;
      for (int n = 0 ; n < nS ; n++) {
         float t = raySphere(V, W, uS[n]);
         if (t > 0. && t < tMin) {
vec3 P=V+t*W;  
            color = shadeSphere(P, uS[n],uSm[n]);
	    tMin = t;

//reflection
 vec3 N = normalize(P - uS[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;
            for (int rn = 0 ; rn < nS ; rn++) {
               float rt = raySphere(P, R, uS[rn]);
               if (rt > 0. && rt < rtMin) {
                  rtMin = rt;
                  rColor = shadeSphere(P + rt * R, uS[rn], uSm[rn]);
               }
            }
            if (rtMin < 10000.)
               color += .3 * rColor;

         }

      }
//end of reflection

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

 let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
   let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   let norm = v => Math.sqrt(dot(v,v));

   let normalize = v => {
      let s = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      return [ v[0]/s, v[1]/s, v[2]/s ];
   }

let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];//closures
let subtract = (a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]];

   let radius = .05;//define size

   let ldData = [ normalize([1,1,1]),
                  normalize([-1,-1,-1]) ];
   S.setUniform('3fv', 'uLd', ldData.flat());
   S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);



   S.setUniform('1f', 'uTime', time);

   // move sphere

   for (let n = 0 ; n < S.nS ; n++) {
 S.sVel[n][0] += 0.001*Math.sin(2*time+n);
S.sVel[n][1] += 0.02*Math.sin(time+n+100);
 S.sVel[n][2] += 0.02*Math.sin(5*time+n);
      for (let i = 0 ; i < 3 ; i++) {
       //  S.sVel[n][i] += .003 * Math.cos(time + (2+i) * n);

         S.sVel[n][i] += .01 * (Math.random() - .5);
         S.sPos[n][i] += .02 * S.sVel[n][i];
      }
      S.sPos[n] = scale(normalize(S.sPos[n]), .7);
   }

   //avoid clashing together

   for (let m = 0 ; m < S.nS ; m++)
   for (let n = 0 ; n < S.nS ; n++)
      if (m != n) {
         let D = subtract(S.sPos[m], S.sPos[n]);
         let d = norm(D);
         if (d < 4 * radius) {
            let t = 4 * radius - d;
            for (let i = 0 ; i < 3 ; i++) {
               S.sPos[m][i] += t * D[i] / d;
               S.sPos[n][i] -= t * D[i] / d;
            }
         }
      }

/*
   let ld0 = normalize([1,1,1]);
   let ld1 = normalize([-1,-1,-1]);

   let ldData = [];
   for (let i = 0 ; i < 3 ; i++)
      ldData.push(ld0[i]);
   for (let i = 0 ; i < 3 ; i++)
      ldData.push(ld1[i]);



   S.setUniform('4fv', 'uS', [ 0,0,0,.5 ]);
   S.setUniform('1f', 'uTime', time);


 for (let n = 0 ; n < S.nS ; n++) {
 for (let i = 0 ; i < 3 ; i++) {
S.sVel[n][i]+=.1*(Math.random()-.5);
S.sPos[n][i]+=.1*S.sVel[n][i];
//S.sPos[n][i]+=.9;
}
S.sPos[n]= normalize(S.sPos[n]);//stop them from moving away


for(let i=0;i<3;i++){
S.sPos[n][i]*=.5;
}

}

*/
   let sData = [];
let smData=[];
   let k = 0;
   for (let n = 0 ; n < S.nS ; n++) {
      //for (let i = 0 ; i < 3 ; i++)
         sData.push(S.sPos[k++]);
      sData.push(.1);//sphere size
	//smData.push([.1,.1,.1,0.,.8,.8,.8,0.,1,1,1,20,0,0,0,0]);
   }

//color
   S.setUniform('3fv', 'uLd', ldData);
   S.setUniform('3fv', 'uLc', [ 1,1,1,lred.value/100,
                                lgreen.value/100,
                                lblue.value/100 ]);
/*   
S.setUniform('3fv', 'uDc', [ red.value/100,
                                green.value/100,
                                blue.value/100 ]);
*/

   S.setUniform('4fv', 'uS', sData.flat());
   S.setUniform('Matrix4fv','uSm',false,S.material.flat());
   //false can be used to transpose matrix

S.setUniform('3fv', 'uBgColor', [ red.value  /100,
                                     green.value/100,
                                     blue.value /100 ]);

   S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
`,
events: `
   ;
`
};

}


