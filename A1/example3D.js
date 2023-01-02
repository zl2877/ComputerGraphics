
rooms.example3D = function() {

lib3D();

description = 'Interactive WebGL<br>on a single square.';

code = {
'explanation': `
   S.html(\`
      Most of the work happens in a fragment shader.
      <p>
      Input to the fragment shader is x,y and time: <code>uPos, uTime</code>
      <p>
      We can also interact by adding information about the cursor: <code>uX,uY</code>
      <p>
      Output at each fragment is: red,green,blue,alpha
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
fragment: `
S.setFragmentShader(\`

   uniform float uTime, uSpace, uX, uY;
   varying vec3 vPos;


vec3 auroralights(float y) {

vec3 aurora =-.1*vec3(.3,3.,5.);      
float s = mix(.6,1., clamp(3.*y-2., 0.,1.));

aurora= mix(aurora, vec3(0.,0.,3.),min(y+1.,1.));

aurora = mix(aurora, vec3(0.1,0.6,0.3), min(y+1.,1.));
     
aurora = mix(aurora, vec3(0.,0.,0.),max(5.*y,0.));


      return aurora;
   }

   void main() {
vec3 p =vPos+vec3(.1*uTime,0.,.1*uTime);
vec3 color=auroralights(p.y+noise(p));//fog like effect
      gl_FragColor = vec4(sqrt(color), 1.);

   }

\`)
`,
render: `
   S.setUniform('1f', 'uTime', time);
   S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
`,
events: `
   onDrag = (x,y) => {
      S.setUniform('1f', 'uX', x);
      S.setUniform('1f', 'uY', y);
   }
   onKeyPress  =k=>S.setUniform('1f','uSpace',k==32);
   onKeyRelease=k=>S.setUniform('1f','uSpace',false);
`
}

}

