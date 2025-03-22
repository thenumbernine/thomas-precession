import {Canvas} from '/js/dom.js';
import {vec3, quat} from '/js/gl-matrix-3.4.1/index.js';
import {getIDs, removeFromParent, show} from '/js/util.js';
import {GLUtil} from '/js/gl-util.js';
import {Mouse3D} from '/js/mouse3d.js';

//TODO add a graph (at least to compare to theoretical values)
//TODO add RK4 integrator

const ids = getIDs();
window.ids = ids;

const urlparams = new URLSearchParams(location.search);

let canvas;
let gl;
let glutil;
let mouse;
let sphereObj;
let axisObj;
let sphereTex;

let axisScale = 1.5;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	glutil.resize();
}

/*
x0 = sinh(gt)/g, x1 = cosh(gt)/g - 1
u0 = cosh(gt), u1=sinh(gt)
a0 = sinh(gt)*g, a1=cosh(gt)*g

reparameterize x1, u0, u1 in terms of x0

cosh^2 - sinh^2 = 1
cosh = sqrt(1 + sinh^2)

(x1*g + 1)^2 - (x0*g)^2 = 1
(x1*g + 1)^2 = (x0*g)^2 + 1
x1*g + 1 = sqrt((x0*g)^2 + 1)
x1*g = sqrt((x0*g)^2 + 1) - 1
x1 = (sqrt((x0*g)^2 + 1) - 1) / g

u0 = sqrt(1 + (x0*g)^2)



x0 = s/sqrt(1-w^2 r^2)
x1 = r/sqrt(1-w^2 r^2) cos(w s)
x2 = r/sqrt(1-w^2 r^2) sin(w s)
x3 = 0

u0 = 1/sqrt(1-w^2 r^2)
u1 = -wr/sqrt(1-w^2 r^2) sin(w s)
u2 = wr/sqrt(1-w^2 r^2) cos(w s)
u3 = 0

a0 = 0
a1 = -w^2 r/sqrt(1-w^2 r^2) cos(w s)
a2 = -w^2 r/sqrt(1-w^2 r^2) sin(w s)
a3 = 0

-Omega = u ^ a = 
[  0							-w^2 r/(1 - w^2 r^2) cos(ws)	-w^2 r/(1 - w^2 r^2) sin(ws)	0 ]
[ w^2 r/(1 - w^2 r^2) cos(ws)			0						w^3 r^2 / (1 - w^2 r^2)			0 ]
[ w^2 r/(1 - w^2 r^2) sin(ws)		-w^3 r^2 / (1 - w^2 r^2)		0							0 ]
[ 0										0							0							0 ]

dS/ds = (-Omega) * S

*/

//story problem:
let omega = .5;
let radius = .5;
let mass = .01;
//electron
//let radius = 2.817940289458e-15;	//m
//let omega = 5151228733459357.28;	// rad/s // angular velocity of the object's orbit
//let mass = 9.10938291e0-31;	//kg
//earth
//let omega = Math.PI * 2 / (60 * 60 * 24);	//rad/s
//let radius = 6371000;	//earth radius in meters
//let mass = 5.972e24;i	//kg

let hbar = 1;		// angular velocity of the rotating object ... how is this different than omega?
let gamma = 1 / Math.sqrt(1 - omega * omega * radius * radius);

let integrators = {
	euler : function(t, x, dt, f) {
		//return x + f(t,x) * t
		//but javascript sucks and doesn't let you overload operators
		return x.add(f(t,x).mul(dt));
	}
};
let integrate = integrators.euler;

function getX(t) {
	return [
		t * gamma, 
		radius * gamma * Math.cos(t * omega), 
		radius * gamma * Math.sin(t * omega), 
		0];
}
function getU(t) {
	return [
		gamma, 
		-omega * radius * gamma * Math.sin(omega * t), 
		omega * radius * gamma * Math.cos(t * omega), 
		0];
}
function getA(t) {
	return [
		0,
		-omega * omega * radius * gamma * Math.cos(t * omega),
		-omega * omega * radius * gamma * Math.sin(t * omega),
		0];
}
function getS(t) {
	return [
		-Math.sqrt(.5) * hbar * omega * radius * gamma * Math.sin(omega * gamma * t),
		Math.sqrt(.5) * hbar * (Math.cos(omega * t) * Math.cos(omega * gamma * t) + gamma * Math.sin(omega * t) * Math.sin(omega * gamma * t)),
		Math.sqrt(.5) * hbar * (Math.sin(omega * t) * Math.cos(omega * gamma * t) - gamma * Math.cos(omega * t) * Math.sin(omega * gamma * t)),
		.5 * hbar];
}

function State(args) {
	this.x = args.x;	//position
	this.u = args.u;	//velocity
	this.S = args.S;	//angular momentum
}
State.prototype = {
	add : function(b) {
		let a = this;
		return new State({
			x : [
				a.x[0] + b.x[0],
				a.x[1] + b.x[1],
				a.x[2] + b.x[2],
				a.x[3] + b.x[3]
			],
			u : [
				a.u[0] + b.u[0],
				a.u[1] + b.u[1],
				a.u[2] + b.u[2],
				a.u[3] + b.u[3]
			],
			S : [
				a.S[0] + b.S[0],
				a.S[1] + b.S[1],
				a.S[2] + b.S[2],
				a.S[3] + b.S[3]
			]
		});
	},
	mul : function(b) {
		let a = this;
		return new State({
			x : [
				a.x[0] * b,
				a.x[1] * b,
				a.x[2] * b,
				a.x[3] * b
			],
			u : [
				a.u[0] * b,
				a.u[1] * b,
				a.u[2] * b,
				a.u[3] * b
			],
			S : [
				a.S[0] * b,
				a.S[1] * b,
				a.S[2] * b,
				a.S[3] * b
			]
		});
	}
};

function getState(t) {
	return new State({
		x : getX(t),
		u : getU(t),
		S : getS(t)
		//omega : [0,0,0,0],
		//theta : [0,0,0,0]
	});
}

let dt = .01; 
let t = 0;
let state = getState(t);
let startTime = undefined;

function update() {

	//update simulation
	//the easy way: implicit
	//if (startTime == undefined) startTime = Date.now();
	//t = (Date.now() - startTime) / 1000;
	//state = getState(t);

	//the hard way: explicit
	state = integrate(t, state, dt, function(t, state) {
		
		// a = du/dtau = -omega * omega * x;
		let negOmegaSq = -omega * omega;
		
		let a = [
			0,
			state.x[1] * negOmegaSq,
			state.x[2] * negOmegaSq,
			state.x[3] * negOmegaSq
		];

		// dS/dtau = (u ^ a) * S
		let dS_dt = [];
		for (let i = 0; i < 4; ++i) {
			dS_dt[i] = 0;
			for (let j = 0; j < 4; ++j) {

				let S_j = j==0 ? -state.S[j] : state.S[j];

				//dS_dt = (u ^ a) * S
				//dS_dt^i = (u ^ a)^ij * S_j
				//dS_dt^i = (u^i a^j - a^j u^i) eta_jk S^k
				dS_dt[i] = dS_dt[i] + (state.u[i] * a[j] - a[i] * state.u[j]) * S_j;
			}
		}

		let dState_dt = new State({
			x : state.u,
			u : a,
			S : dS_dt
			//TODO ???
			//omega : [omega, 0, 0, 0],
			//theta : [0, 0, 0, 0]
		});
		return dState_dt;
	});
	t += dt;

	/*
	J_uv = U^a S^b e_abuv
	so J_uv is the tangent space to the rotation axis ... in the spatial space
	*/
	
	let omegaAxis = vec3.fromValues(
		state.S[1] * 5 / (2 * mass * radius * radius),
		state.S[2] * 5 / (2 * mass * radius * radius),
		state.S[3] * 5 / (2 * mass * radius * radius));
	//and with the angular velocity, integrate by time (what kind of time?) to get the angle itself
	let omegaAngle = vec3.length(omegaAxis);
	if (Math.abs(omegaAngle) > 1e-20) {
		let omegaQ = quat.create();
		vec3.scale(omegaAxis, omegaAxis, 1/omegaAngle);
		omegaAngle *= dt;
		quat.setAxisAngle(omegaQ, omegaAxis, omegaAngle);
		quat.multiply(sphereObj.angle, sphereObj.angle, omegaQ);
		quat.copy(axisObj.angle, sphereObj.angle);

		axisObj.attrs.vertex.buffer.data[0] = axisScale*omegaAxis[0];
		axisObj.attrs.vertex.buffer.data[1] = axisScale*omegaAxis[1];
		axisObj.attrs.vertex.buffer.data[2] = axisScale*omegaAxis[2];
		axisObj.attrs.vertex.buffer.data[3] = -axisScale*omegaAxis[0];
		axisObj.attrs.vertex.buffer.data[4] = -axisScale*omegaAxis[1];
		axisObj.attrs.vertex.buffer.data[5] = -axisScale*omegaAxis[2];
		axisObj.attrs.vertex.buffer.updateData();
	}
//console.log(state.S[0], state.S[1], state.S[2], state.S[3]);	
	glutil.draw();
	requestAnimationFrame(update);
}

canvas = Canvas({
	style : {
		left : 0,
		top : 0,
		position : 'absolute',
		userSelect : 'none',
	},
	prependTo : document.body,
});

try {
	glutil = new GLUtil({canvas:canvas});
	gl = glutil.context;
} catch (e) {
	panel.remove();
	removeFromParent(canvas);
	show(ids.webglfail);
	throw e;
}
gl.enable(gl.DEPTH_TEST);
glutil.view.pos[0] = -6;
glutil.view.fovY = 45;

//x fwd, z up, y left:
let SQRT_1_2 = Math.sqrt(1/2);
quat.mul(glutil.view.angle, /*90' x*/[SQRT_1_2,0,0,SQRT_1_2], /*90' -y*/[0,-SQRT_1_2,0,SQRT_1_2]);

sphereTex = new glutil.Texture2D({
	flipY : true,
	generateMipmap : true,
	magFilter : gl.LINEAR,
	minFilter : gl.LINEAR_MIPMAP_LINEAR,
	url : 'BlueMarble.png',
	onload : () => {
		let sphereShader = new glutil.Program({
			vertexCode : `
in vec3 vertex;
in vec3 tc;
out vec2 vtc;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	vec4 eyePos = mvMat * vec4(vertex, 1.);
	vtc = tc.xy;
	gl_Position = projMat * eyePos;
}
`,
			fragmentCode : `
in vec2 vtc;
uniform sampler2D tex;
out vec4 fragColor;
void main() {
	fragColor = texture(tex, vtc);
}
`,
		});

		let w = 20;
		let h = 10;
		let sphereVtxArray = new Float32Array(3 * 2 * 3 * w * h);
		let sphereTCArray = new Float32Array(3 * 2 * 3 * w * h);
		let radius = 1;
		let writePt = function(ofs, u, v) {
			let phi = 2 * Math.PI * u;
			let theta = Math.PI * v;
			//should be x y z but im too lazy to correct the camera for y being up
			let z = radius * Math.cos(phi) * Math.sin(theta);
			let x = radius * Math.sin(phi) * Math.sin(theta);
			let y = radius * Math.cos(theta);
			sphereVtxArray[0 + ofs] = x;
			sphereVtxArray[1 + ofs] = y;
			sphereVtxArray[2 + ofs] = z;
			sphereTCArray[0 + ofs] = u;
			sphereTCArray[1 + ofs] = 1 - v;
			sphereTCArray[2 + ofs] = 0;
		};
		for (let v = 0; v < h; v++) {
			for (let u = 0; u < w; u++) {
				writePt(0 + 18 * (u + w * v), u/w, v/h);
				writePt(3 + 18 * (u + w * v), (u+1)/w, v/h);
				writePt(6 + 18 * (u + w * v), (u+1)/w, (v+1)/h);
				writePt(9 + 18 * (u + w * v), (u+1)/w, (v+1)/h);
				writePt(12 + 18 * (u + w * v), u/w, (v+1)/h);
				writePt(15 + 18 * (u + w * v), u/w, v/h);
			}
		}
		
		sphereObj = new glutil.SceneObject({
			mode : gl.TRIANGLES,
			shader : sphereShader,
			static : false,
			uniforms : { 
				tex : 0
			},
			attrs : { 
				vertex : new glutil.ArrayBuffer({
					data : sphereVtxArray
				}),
				tc : new glutil.ArrayBuffer({
					data : sphereTCArray
				})
			},
			texs : [sphereTex]
		});

		axisObj = new glutil.SceneObject({
			mode : gl.LINES,
			shader : new glutil.Program({
				vertexCode : `
in vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex, 1.);
}
`,
				fragmentCode : `
out vec4 fragColor;
void main() {
	fragColor = vec4(0., 1., 1., 1.); 
}
`,
			}),
			attrs : {
				vertex : new glutil.ArrayBuffer({
					data : [0,0,axisScale*radius, 0,0,-axisScale*radius],
					usage : gl.DYNAMIC_DRAW
				})
			},
			static : false
		});

		let lastMouseRot = quat.create();	
		let tmpQ = quat.create();
		let tmpV = vec3.create();
		mouse = new Mouse3D({
			pressObj : canvas,
			move : function(dx,dy) {
				let rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
				quat.setAxisAngle(tmpQ, [-dy, -dx, 0], rotAngle);
				//mat4.translate(glutil.scene.mvMat, glutil.scene.mvMat, [10*dx/canvas.width, -10*dy/canvas.height, 0]);

				//put tmpQ into the frame of glutil.view.angle, so we can rotate the view vector by it
				//  lastMouseRot = glutil.view.angle-1 * tmpQ * glutil.view.angle
				// now newViewAngle = glutil.view.angle * tmpQ = lastMouseRot * glutil.view.angle
				// therefore lastMouseRot is the global transform equivalent of the local transform of tmpQ
				quat.mul(lastMouseRot, glutil.view.angle, tmpQ);
				quat.conjugate(tmpQ, glutil.view.angle);
				quat.mul(lastMouseRot, lastMouseRot, tmpQ);

				vec3.copy(tmpV, glutil.view.pos);
				let posDist = vec3.length(tmpV);
				vec3.transformQuat(tmpV, tmpV, lastMouseRot);
				vec3.normalize(tmpV, tmpV);
				vec3.scale(tmpV, tmpV, posDist);
				vec3.copy(glutil.view.pos, tmpV);
				
				//RHS apply so it is relative to current view 
				//newViewAngle := glutil.view.angle * tmpQ
				quat.mul(glutil.view.angle, lastMouseRot, glutil.view.angle);
				quat.normalize(glutil.view.angle, glutil.view.angle);

			},
			zoom : function(dz) {
				glutil.view.fovY *= Math.exp(-.0003 * dz);
				glutil.view.fovY = Math.clamp(glutil.view.fovY, 1, 179);
				glutil.updateProjection();
			}
		});
		resize();
		update();
	},
});
