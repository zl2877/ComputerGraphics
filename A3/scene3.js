
rooms.scene3 = function () {

   lib3D2();

   description = `<b>Scene 3</b>
   <p>For this scene I created many vases using revolution, and a ribbon at the back using extrusions
   <p>The vase shape can be changed by slider
   <p><b>Change color</b> </p>
                  <br> <input type=range id=red   value= 20> red
                  <br> <input type=range id=green value=60> green
                  <br> <input type=range id=blue  value=100> blue
                  <p> <input type=range id=bottleShape value=1  min="1" max="3">Vase shape
              `;

   code = {
      'init': `


























   S.material = [
      [.5,.1,.1,0,  .5,.2,.2,0,  1,1,1,5,    0,0,0,0], //S.redPlastic 
      [.1,.3,.1,0,  .2,.5,.2,0,   2,2,2,20,   0,0,0,0], //S.greenPlastic 
      [.1,.2,.5,0,  .2,.2,.5,0,  2,2,2,20, 0,0,0,0], //S.bluePlastic  
      [.2,.2,.2,0,  .5,.5,.5,0, 1,.6,.1,6,  0,0,0,0], //S.whitePlastic  
      [.25,.05,.25,0,  .5,.5,.5,0, 1,1,1,5,    0,0,0,0], //S.purplePlastic  
      /*
      [.5,.1,.1,0,  .5,.2,.2,0,  2,2,2,20,  0,0,0,0];//S.redPlastic   
      [.1,.3,.1,0,  .2,.5,.2,0,  2,2,2,20,  0,0,0,0];//S.greenPlastic  
       [.1,.2,.5,0,  .2,.2,.9,0,  2,2,2,20,  0,0,0,0];//S.bluePlastic  
      [.2,.2,.2,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];//S.whitePlastic  
      [.25,.05,.25,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];//S.purplePlastic  
      
      [.1,.1,.1,0,     .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // SILVER
      [.25,0,0,0,      .5,0,0,0,    2,2,2,20,   0,0,0,0], // PLASTIC
      [.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,3, 0,0,0,0], // COPPER
      [.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0], // GOLD
      [.05,.05,.05,0,  .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // LEAD
      */
   ];
   
   S.nM = S.material.length;

   // A SQUARE IS A TRIANGLE MESH WITH JUST TWO TRIANGLES

   S.squareMesh = [ -1, 1, 0,  0,0,1,  0,1,
                     1, 1, 0,  0,0,1,  1,1,
                    -1,-1, 0,  0,0,1,  0,0,
                     1,-1, 0,  0,0,1,  1,0 ];

   // GLUE TOGETHER TWO MESHES TO CREATE A SINGLE MESH

   let glueMeshes = (a,b) => {
      let mesh = a.slice();
      mesh.push(a.slice(a.length - S.VERTEX_SIZE, a.length));
      mesh.push(b.slice(0, S.VERTEX_SIZE));
      mesh.push(b);
      return mesh.flat();
   }

   let add      = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]; //to add
   let subtract = (a,b) => [ b[0]-a[0], b[1]-a[1], b[2]-a[2] ];//to subtract 
   let cross    = (a,b) => [ a[1] * b[2] - a[2] * b[1],
                             a[2] * b[0] - a[0] * b[2],
                             a[0] * b[1] - a[1] * b[0] ];
   let norm = a => Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
   let normalize = a => {
      let s = norm(a);
      return s < .00001 ? [0,0,0] : [ a[0] / s, a[1] / s, a[2] / s ];
   }

   // GIVEN A FUNCTION THAT MAPS (u,v) TO point AND normal,
   // AND GIVEN A MESH RESOLUTION, CREATE A PARAMETRIC MESH

   let uvMesh = (f, nu,nv, data) => {
      let mesh = [];

/////// YOU NEED TO IMPLEMENT THE FOLLOWING SECTIONS ///////////////////////////////////////////////////////////////////////////////////

      // CREATE AN ARRAY OF nu+1 X nv+1 VERTICES
/*
           v---v---v
           |   |   |
           v---v---v
           |   |   |
           v---v---v
           |   |   |
           v---v---v
*/

// CREATE AN ARRAY OF nu X nv FACE NORMALS
/*
           D---C---v
           | f |   |
           A---B---v  f = (B-A) X (C-B) +
           |   |   |      (C-B) X (D-C) +
           v---v---v      (D-C) X (A-D) +
           |   |   |      (A-D) X (B-A)
           v---v---v
*/

let vertices = [];
let faceNormal = [];
let v=0;
let u=0;

for (v=0; v <= nv; v++) {
   for (u = 0; u <= nu; u++) {
      vertices.push(f(u/nu, v/nv, data).slice(0, 3));
   }
}

for (v = 0; v < nv; v++) {
   for (u = 0; u < nu; u++) {
      let A = vertices[u+v*(nu+1)];
      let B = vertices[u+1+v*(nu+1)];
      let C = vertices[u +1+(nu + 1)*(v+1)];
      let D = vertices[u+(v+1) * (nu + 1)];
      let f = add(cross((subtract(A, B)), subtract(B, C)), //(B-A) X (C-B)
              add(cross(subtract(B, C), subtract(C, D)), //(C-B) X (D-C)
              add(cross(subtract(C, D), subtract(D, A)), //(D-C) X (A-D)
              cross(subtract(D, A), subtract(A, B))))); //(A-D) X (B-A)
      faceNormal.push(f);
   }
}

      // SUM THE 4 ADJOINING FACE NORMALS TO COMPUTE EACH VERTEX NORMAL
/*
           d---c---v
           |f2 |f3 |
           a---N---v   N = normalize(f0 + f1 + f2 + f3)
           |f0 |f1 |
           v---v---v
           | f | f |
           v---v---v
*/

let vertexNormal = [];
let N;
let f0,f1,f2,f3;
let f0_u, f0_v;
let f1_u, f1_v;
let f2_u, f2_v;
let f3_u, f3_v;

		for (v = 0; v <= nv; v++) {
			for (u = 0; u <= nu; u++) {
				f0_u = u - 1;
            f0_v = v - 1;

				f1_u = u;
            f1_v = v - 1;

				f2_v = v;
				f2_u = u - 1;

				f3_u = u;
				f3_v = v;
				
				if (u == 0) {
					f0_u = nu - 1;
					f2_u = nu - 1;
				}
				if (u == nu) {
					f1_u = 0;
					f3_u = 0;
				}
				if (v == 0) {
					f0_u = nu - 1 - f3_u;
					f1_u = nu - 1 - f2_u;
					f0_v = 0;
					f1_v = 0;
				}
				if (v == nv) {
					f2_u = nu - 1 - f1_u;
					f3_u = nu - 1 - f0_u;
					f2_v = v - 1;
					f3_v = v - 1;
				}
				f0 = faceNormal[f0_v * nu + f0_u];
				f1 = faceNormal[f1_v * nu + f1_u];
				f2 = faceNormal[f2_v * nu + f2_u];
				f3 = faceNormal[f3_v * nu + f3_u];

            //N = normalize(f0 + f1 + f2 + f3)
            N=normalize(add((add(f0,f1)),add(f2,f3)));
				vertexNormal.push(N);
			}
		}

      // BUILD THE MESH BY GLUEING TOGETHER ROWS OF TRIANGLE STRIPS
/*
        Don't try to build a flat array here.
        Make this an array of arrays, where each vertex is its own array.
        In particular, use mesh.push() rather than mesh.concat().
*/

for (v= 0; v < nv; v++) {
   mesh.push(vertices[v * (nu + 1)]);
   mesh.push(vertexNormal[v * (nu + 1)]);
   mesh.push([0, v]);
   for (u = 0; u <= nu; u++) {

      /*
      let A = vertices[u+v*(nu+1)];
      let B = vertices[u+1+v*(nu+1)];
      let C = vertices[u +1+(nu + 1)*(v+1)];
      let D = vertices[u+(v+1) * (nu + 1)];
      */
      mesh.push(vertices[u+(v+1)*(nu+1)]);
      mesh.push(vertexNormal[u+(v+1)*(nu+1)]);
      mesh.push([u,v+1]);

      mesh.push(vertices[u+v*(nu+1)]);
      mesh.push(vertexNormal[u+v*(nu+1)]);
      mesh.push([u, v]);
   }
   mesh.push(vertices[nu+(v+1)*(nu+1)]);
   mesh.push(vertexNormal[nu+(v+1)*(nu+1)]);
   mesh.push([nu, v+1]);
}
      // RETURN THE FLATTENED ARRAY
/*
        Finally, just flatten everything using the .flat() method.
*/
        return mesh.flat();
   }

   S.uvMesh = uvMesh;

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // CREATE A UNIT SPHERE PARAMETRIC MESH

   S.sphereMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = Math.PI * v - Math.PI/2;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      let cv = Math.cos(phi);
      let sv = Math.sin(phi);
      return [cu * cv, su * cv, sv,
              cu * cv, su * cv, sv,
              u, v];
   }, 20, 10);

   // CREATE A UNIT TORUS PARAMETRIC MESH

   S.torusMesh = uvMesh((u,v,r) => {
      let theta = 2 * Math.PI * u;
      let phi   = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      let cv = Math.cos(phi);
      let sv = Math.sin(phi);
      return [cu * (1 + r * cv), su * (1 + r * cv), r * sv,
              cu * cv, su * cv, sv,
              u, v];
   }, 20, 10, .4);

   // CREATE A UNIT DISK PARAMETRIC MESH

   S.diskMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi   = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [v * cu, v * su, 0,  0, 0, 1,   u, v];
   }, 20, 2);

   // CREATE A UNIT OPEN TUBE PARAMETRIC MESH

   S.tubeMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi   = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [cu, su, 2 * v - 1,   cu, su, 0,   u, v];
   }, 20, 2);

   // TRANSFORM A MESH BY A MATRIX ON THE CPU

   let transformMesh = (mesh, matrix) => {
      let result = [];
      let IMT = matrixTranspose(matrixInverse(matrix));
      for (let n = 0 ; n < mesh.length ; n += S.VERTEX_SIZE) {
         let V = mesh.slice(n, n + S.VERTEX_SIZE);
         let P  = V.slice(0, 3);
         let N  = V.slice(3, 6);
         let UV = V.slice(6, 8);
         P = matrixTransform(matrix, [P[0], P[1], P[2], 1]);
         N = matrixTransform(IMT,    [N[0], N[1], N[2], 0]);
         result.push(P[0],P[1],P[2], N[0],N[1],N[2], UV);
      }
      return result.flat();
   }

   // A CYLINDER MESH IS A TUBE WITH TWO DISK END-CAPS GLUED TOGETHER

   let end0 = transformMesh(S.diskMesh, matrixTranslate([0,0,1]));
   let end1 = transformMesh(end0      , matrixRotx(Math.PI));
   S.cylinderMesh = glueMeshes(S.tubeMesh, glueMeshes(end0, end1));

   // A CUBE MESH IS SIX TRANSFORMED SQUARE MESHES GLUED TOGETHER

   let face0 = transformMesh(S.squareMesh, matrixTranslate([0,0,1]));
   let face1 = transformMesh(face0,        matrixRotx( Math.PI/2));
   let face2 = transformMesh(face0,        matrixRotx( Math.PI  ));
   let face3 = transformMesh(face0,        matrixRotx(-Math.PI/2));
   let face4 = transformMesh(face0,        matrixRoty(-Math.PI/2));
   let face5 = transformMesh(face0,        matrixRoty( Math.PI/2));
   S.cubeMesh = glueMeshes(face0,
                glueMeshes(face1,
                glueMeshes(face2,
                glueMeshes(face3,
                glueMeshes(face4,
                           face5)))));

   // DRAW A SINGLE MESH.

   S.drawMesh = (mesh, matrix, materialIndex) => {
      let gl = S.gl;
      if (! S.gl.bufferData)
         return;
      S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
      S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
      S.setUniform('Matrix4fv', 'uMaterial', false, S.material[materialIndex]);
      S.gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      S.gl.drawArrays(mesh.isTriangles ? S.gl.TRIANGLES
                                       : S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
   }

   let evalCubicSpline = (splineMatrix, P, t) => {
      let splineValue = P => {
         let C = matrixTransform(splineMatrix, P);
         return t*t*t * C[0] + t*t * C[1] + t * C[2] + C[3];
      }

      // THE VALUE AT A KEY CAN BE EITHER A NUMBER OR AN OBJECT

      if (Number.isFinite(P[0]))    // SPECIAL CASE: THE VALUE
         return splineValue(P);     // AT THE KEY IS A NUMBER.

      let value = {};
      for (let k in P[0])
         value[k] = splineValue([ P[0][k], P[1][k], P[2][k], P[3][k] ]);
      return value;
   }

   let CatmullRomMatrix = [
     -1/2,  1  , -1/2, 0,
      3/2, -5/2,  0  , 1,
     -3/2,  2  ,  1/2, 0,
      1/2, -1/2,  0  , 0,
   ];

   S.CatmullRomFunction = (keys, n, t) => {
      let mm = n => Math.max(0, Math.min(keys.length - 1, n));
      let a = keys[mm(n-1)];
      let b = keys[mm(n  )];
      let c = keys[mm(n+1)];
      let d = keys[mm(n+2)];
      return evalCubicSpline(CatmullRomMatrix, [a,b,c,d], t);
   }

   S.evalSpline = (keys, f, splineFunction) => {
      let T = Math.max(0, Math.min(.9999, f)) * (keys.length - 1);
      return splineFunction(keys, T >> 0, T % 1);
   }
 
   // CREATE A SURFACE OF REVOLUTION MESH

   S.createRevolutionMesh = (nu,nv,keys) => S.uvMesh((u,v,keys) => {
      let theta = 2 * Math.PI * u;
      let cos = Math.cos(theta);
      let sin = Math.sin(theta);

      let zr  = S.evalSpline(keys, v, S.CatmullRomFunction);

      return [
         zr.r * cos, zr.r * sin, zr.z,
         0,0,0,                // NORMAL WILL BE COMPUTED LATER IN uvMesh().
         u, v
      ];
   }, nu, nv, keys);

   S.createExtrusionMesh = (nu,nv,data) => {

      let radius   = data.radius;
      let profile  = data.profile;
      let path     = data.path;
      let profileSpline = u => S.evalSpline(profile, u, S.CatmullRomFunction);
      let pathSpline    = v => S.evalSpline(path   , v, S.CatmullRomFunction);

      let m = new Matrix(),
          p = pathSpline(0),
          q = pathSpline(0.001);
/*


      /////// YOU NEED TO IMPLEMENT THE FOLLOWING SECTION ////////////////////////////////////////////////////////

      Z = NORMALIZE(q - p)
      X = A VECTOR NOT ALIGNED WITH Z

      // TO FIND A REASONABLE INITIAL VALUE FOR X:

         xx = Z[0]*Z[0]
         yy = Z[1]*Z[1]
         zz = Z[2]*Z[2]

         if xx < yy && xx < zz then X = [1,0,0]
         if yy < xx && yy < zz then X = [0,1,0]
         if zz < xx && zz < yy then X = [0,0,1]
*/

let Z,X;
let xx,yy,zz;
//Z = NORMALIZE(q - p)
//Z = normalize(subtract(q-p));
Z = normalize([q.x - p.x, q.y - p.y, q.z - p.z]);
xx = Z[0]*Z[0];
yy = Z[1]*Z[1];
zz = Z[2]*Z[2];
if (xx < yy && xx < zz){
   X = [1,0,0];
}
if (yy < xx && yy < zz){
   X = [0,1,0];
}
if (zz < xx && zz < yy){ 
   X = [0,0,1];
}
return S.uvMesh((u,v) => {
   p = pathSpline(v - .001);
   q = pathSpline(v + .001);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
         /////// YOU NEED TO IMPLEMENT THE FOLLOWING SECTION ///////

         Z = NORMALIZE(q - p)
         Y = NORMALIZE( CROSS (Z, X) )
         X = NORMALIZE( CROSS (Y, Z) )
         m = X Y Z p
*/

let Z,Y;
Z =normalize([q.x - p.x, q.y - p.y, q.z - p.z]);
Y = normalize(cross(Z, X));
X = normalize(cross(Y, Z));
m.set([X,Y,Z,p.x, p.y, p.z].flat());
m.set([X, 0, Y, 0, Z, 0, p.x, p.y, p.z, 0].flat());
//m.flat();


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
   p = profileSpline(u);
   let P = m.transform([ radius * p.x, radius * p.y, radius * p.z ]);
   return [
      P[0],P[1],P[2],
      0,0,0,            // NORMAL WILL BE COMPUTED LATER IN uvMesh().
      u,v
   ];

}, nu, nv);
};
`,

      fragment: `
S.setFragmentShader(\`
   const int nL = \` + S.nL + \`;
   const int nM = \` + S.nM + \`;
   uniform vec3 uBgColor;
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];
   uniform mat4 uMaterial;
   varying vec3 vPos, vNor;

   void main() {
      vec3 N = normalize(vNor);
      vec3  ambient  = uMaterial[0].rgb;
      vec3  diffuse  = uMaterial[1].rgb;
      vec3  specular = uMaterial[2].rgb;
      float p        = uMaterial[2].a;
      vec3 c = mix(ambient, uBgColor, .3);
      for (int l = 0 ; l < nL ; l++) {
         vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
         c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                      + specular * pow(max(0., R.z), p));
      }
      gl_FragColor = vec4(c, 1.);
   }
\`);
`,
      vertex: `
S.setVertexShader(\`
   attribute vec3 aPos, aNor;
   uniform   mat4 uMatrix, uInvMatrix, uProject;
   varying   vec3 vPos, vNor;

   void main() {
      vPos = (uProject * uMatrix * vec4(aPos, 1.)).xyz;
      vNor = (vec4(aNor, 0.) * uInvMatrix).xyz;
      gl_Position = vec4(vPos.xy, -.01 * vPos.z, 1.);
   }
\`)
`,
      render: `


   S.revolutionMesh = S.createRevolutionMesh(20, 32, [
      //{z:-2 , r:.5 },
      {z:-1.1 , r:0 },
      {z:-1.05  , r:0.3  },
      //{z:-.1, r:.1 },
      {z:-.9 , r:.2/bottleShape.value  },
      {z:-.5 , r:.5 },
      {z: .3 , r:.1 },
      {z: .9 , r:.18/bottleShape.value  },
      
      {z: 1  , r:0.19  },
      
   ]);

   S.revolutionMesh2 = S.createRevolutionMesh(20, 32, [
      {z:-1.1 , r:0 },
      {z:-1.07  , r:0.2  },
      
      {z:-.95 , r:.2 },
      {z:-.7 , r:.5 },
      {z:-.2, r:.7/bottleShape.value  },
      {z: .3 , r:.5 },
      {z: .7 , r:.1 },
      {z: 1  , r:0.1  },
      
   ]);

   S.revolutionMesh3 = S.createRevolutionMesh(20, 32, [
      {z:-1.1 , r:0 },
      {z:-1.1  , r:.4  },
      
      {z:-.8 , r:.15/bottleShape.value  },
      {z:-.5 , r:.4  },
      {z:-.2, r:.15/bottleShape.value  },
      {z: .1 , r:.4 },
      {z: .4 , r:.15/bottleShape.value  },
      {z: .7 , r:.4 },
      {z: 1.0 , r:.1 },
      {z: 1.1 , r:0.1  },
      
   ]);


   S.revolutionMesh4 = S.createRevolutionMesh(20, 32, [
      {z:-1.1 , r:0 },
      {z:-1.1  , r:.3  },
   

      {z: .7 , r:.4/bottleShape.value },
      {z: 1.0 , r:.2 },
      {z: 1.1 , r:0.2  },
      
   ]);

   let extrusionData = {
      radius: 0.25,
      profile: [
       //  {x:-1 , y:-1 , z: 1},
       //  {x: 1 , y:-1 , z: 1},
         {x: 1 , y: 1 , z: 1},
         {x:-1 , y: 1 , z: 1},
        
         
         
      ],
      path: [
         {x:-5, y:-2, z: -1},
         {x: -2, y:-1, z: 0},
         {x: 1, y: 1, z: 0},
         {x:4, y: 3, z: 2},

         
      ]
   };

   S.extrusionMesh = S.createExtrusionMesh(38, 28, extrusionData);
   
   

   // SET THE PROJECTION MATRIX BASED ON CAMERA FOCAL LENGTH

   let fl = .0;
   S.setUniform('Matrix4fv', 'uProject', false,
      [1,0,0,0, 0,1,0,0, 0,0,1,-1/fl, 0,0,0,1]);

   // SPECIFY SCENE LIGHTING

   S.nL = 2;
   S.setUniform('3fv', 'uLd', [ .57,.57,.57, -.57,-.57,-.57 ]);
   S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);
   //S.setUniform('3fv', 'uBgColor', [ .89,.81,.75 ]);
   S.setUniform('3fv', 'uBgColor', [ red.value  /100,
                                     green.value/100,
                                     blue.value /100 ]);
   // RENDER THE SCENE

   let m = new Matrix();
   m.save();
      m.scale(.3);
      m.translate(1.5, -1.5, .07);
      m.rotx(Math.PI/8 * Math.sin(time)+5);
      m.roty(Math.PI/8 * Math.cos(time));
      S.drawMesh(S.revolutionMesh3, m.get(), 4);
     // S.drawMesh(S.extrusionMesh, m.get(), 1);
   m.restore();

   m.save();
   m.scale(.3);
   m.translate(-1.5, -1.5, .07);
   m.rotx(Math.PI/8 * Math.sin(time)+5);
   m.roty(Math.PI/8 * Math.cos(time));
   S.drawMesh(S.revolutionMesh2, m.get(), 0);
  // S.drawMesh(S.extrusionMesh, m.get(), 1);
m.restore();

m.save();
m.scale(.3);
m.translate(1.5, 1.5, .07);
m.rotx(Math.PI/8 * Math.sin(time)+5);
m.roty(Math.PI/8 * Math.cos(time));
S.drawMesh(S.revolutionMesh, m.get(), 1);
// S.drawMesh(S.extrusionMesh, m.get(), 1);
m.restore();

m.save();
m.scale(.3);
m.translate(-1.5, 1.5, .07);
m.rotx(Math.PI/8 * Math.sin(time)+5);
m.roty(Math.PI/8 * Math.cos(time));
S.drawMesh(S.revolutionMesh4, m.get(), 2);
// S.drawMesh(S.extrusionMesh, m.get(), 1);
m.restore();

   m.save();
   m.scale(.3);
   m.translate(.0, .0, -2.);
   m.rotx(Math.PI/8 * Math.sin(time));
   m.roty(Math.PI/18 * Math.cos(time));
 //  S.drawMesh(S.revolutionMesh, m.get(), 0);
   S.drawMesh(S.extrusionMesh, m.get(), 3);
   m.restore();
`,
      events: `
   ;
`
   };

}

