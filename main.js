//TODO add a graph (at least to compare to theoretical values)
//TODO add RK4 integrator

var canvas;
var gl;
var renderer;
var mouse;
var sphereObj;
var axisObj;
var axisScale = 1.5;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	renderer.resize();
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
var omega = .5;
var radius = .5;
var mass = .01;
//electron
//var radius = 2.817940289458e-15;	//m
//var omega = 5151228733459357.28;	// rad/s // angular velocity of the object's orbit
//var mass = 9.10938291e0-31;	//kg
//earth
//var omega = Math.PI * 2 / (60 * 60 * 24);	//rad/s
//var radius = 6371000;	//earth radius in meters
//var mass = 5.972e24;i	//kg

var hbar = 1;		// angular velocity of the rotating object ... how is this different than omega?
var gamma = 1 / Math.sqrt(1 - omega * omega * radius * radius);

var integrators = {
	euler : function(t, x, dt, f) {
		//return x + f(t,x) * t
		//but javascript sucks and doesn't let you overload operators
		return x.add(f(t,x).mul(dt));
	}
};
var integrate = integrators.euler;

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
		var a = this;
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
		var a = this;
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

var dt = .01; 
var t = 0;
var state = getState(t);
var startTime = undefined;

function update() {

	//update simulation
	//the easy way: implicit
	//if (startTime == undefined) startTime = Date.now();
	//t = (Date.now() - startTime) / 1000;
	//state = getState(t);

	//the hard way: explicit
	state = integrate(t, state, dt, function(t, state) {
		
		// a = du/dtau = -omega * omega * x;
		var negOmegaSq = -omega * omega;
		
		var a = [
			0,
			state.x[1] * negOmegaSq,
			state.x[2] * negOmegaSq,
			state.x[3] * negOmegaSq
		];

		// dS/dtau = (u ^ a) * S
		var dS_dt = [];
		for (var i = 0; i < 4; ++i) {
			dS_dt[i] = 0;
			for (var j = 0; j < 4; ++j) {

				var S_j = j==0 ? -state.S[j] : state.S[j];

				//dS_dt = (u ^ a) * S
				//dS_dt^i = (u ^ a)^ij * S_j
				//dS_dt^i = (u^i a^j - a^j u^i) eta_jk S^k
				dS_dt[i] = dS_dt[i] + (state.u[i] * a[j] - a[i] * state.u[j]) * S_j;
			}
		}

		var dState_dt = new State({
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
	
	var omegaAxis = vec3.fromValues(
		state.S[1] * 5 / (2 * mass * radius * radius),
		state.S[2] * 5 / (2 * mass * radius * radius),
		state.S[3] * 5 / (2 * mass * radius * radius));
	//and with the angular velocity, integrate by time (what kind of time?) to get the angle itself
	var omegaAngle = vec3.length(omegaAxis);
	if (Math.abs(omegaAngle) > 1e-20) {
		var omegaQ = quat.create();
		vec3.scale(omegaAxis, omegaAxis, 1/omegaAngle);
		omegaAngle *= dt;
		quat.setAxisAngle(omegaQ, omegaAxis, omegaAngle);
		quat.multiply(sphereObj.angle, sphereObj.angle, omegaQ);
		quat.copy(axisObj.angle, sphereObj.angle);

		axisObj.attrs.vertex.data[0] = axisScale*omegaAxis[0];
		axisObj.attrs.vertex.data[1] = axisScale*omegaAxis[1];
		axisObj.attrs.vertex.data[2] = axisScale*omegaAxis[2];
		axisObj.attrs.vertex.data[3] = -axisScale*omegaAxis[0];
		axisObj.attrs.vertex.data[4] = -axisScale*omegaAxis[1];
		axisObj.attrs.vertex.data[5] = -axisScale*omegaAxis[2];
		axisObj.attrs.vertex.updateData();
	}
//console.log(state.S[0], state.S[1], state.S[2], state.S[3]);	
	renderer.draw();
	requestAnimFrame(update);
}

$(document).ready(init1);

function init1() {
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection();

	try {
		renderer = new GL.CanvasRenderer({canvas:canvas});
		gl = renderer.context;
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	gl.enable(gl.DEPTH_TEST);
	renderer.view.pos[0] = -6;
	renderer.view.fovY = 45;
	
	//x fwd, z up, y left:
	var SQRT_1_2 = Math.sqrt(1/2);
	quat.mul(renderer.view.angle, /*90' x*/[SQRT_1_2,0,0,SQRT_1_2], /*90' -y*/[0,-SQRT_1_2,0,SQRT_1_2]);

	sphereTex = new GL.Texture2D({
		context : gl,
		flipY : true,
		generateMipmap : true,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR_MIPMAP_LINEAR,
		url : 'BlueMarble.png',
		onload : function() {
			init2(this);
		}
	});
}

function init2(sphereTex) {
	var sphereShader = new GL.ShaderProgram({
		context : gl,
		vertexCodeID : 'sphere-vsh',
		fragmentCodeID : 'sphere-fsh'
	});

	var w = 20;
	var h = 10;
	var sphereVtxArray = new Float32Array(3 * 2 * 3 * w * h);
	var sphereTCArray = new Float32Array(3 * 2 * 3 * w * h);
	var radius = 1;
	var writePt = function(ofs, u, v) {
		var phi = 2 * Math.PI * u;
		var theta = Math.PI * v;
		//should be x y z but im too lazy to correct the camera for y being up
		var z = radius * Math.cos(phi) * Math.sin(theta);
		var x = radius * Math.sin(phi) * Math.sin(theta);
		var y = radius * Math.cos(theta);
		sphereVtxArray[0 + ofs] = x;
		sphereVtxArray[1 + ofs] = y;
		sphereVtxArray[2 + ofs] = z;
		sphereTCArray[0 + ofs] = u;
		sphereTCArray[1 + ofs] = 1 - v;
		sphereTCArray[2 + ofs] = 0;
	};
	for (var v = 0; v < h; v++) {
		for (var u = 0; u < w; u++) {
			writePt(0 + 18 * (u + w * v), u/w, v/h);
			writePt(3 + 18 * (u + w * v), (u+1)/w, v/h);
			writePt(6 + 18 * (u + w * v), (u+1)/w, (v+1)/h);
			writePt(9 + 18 * (u + w * v), (u+1)/w, (v+1)/h);
			writePt(12 + 18 * (u + w * v), u/w, (v+1)/h);
			writePt(15 + 18 * (u + w * v), u/w, v/h);
		}
	}
	
	sphereObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		mode : gl.TRIANGLES,
		shader : sphereShader,
		static : false,
		uniforms : { 
			tex : 0
		},
		attrs : { 
			vertex : new GL.ArrayBuffer({
				context : gl,
				data : sphereVtxArray
			}),
			tc : new GL.ArrayBuffer({
				context : gl,
				data : sphereTCArray
			})
		},
		texs : [sphereTex]
	});

	axisObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		mode : gl.LINES,
		shader : new GL.ShaderProgram({
			context : gl,
			vertexCodeID : 'axis-vsh',
			fragmentCodeID : 'axis-fsh'
		}),
		attrs : {
			vertex : new GL.ArrayBuffer({
				context : gl,
				data : [0,0,axisScale*radius, 0,0,-axisScale*radius],
				usage : gl.DYNAMIC_DRAW
			})
		},
		static : false
	});

	var lastMouseRot = quat.create();	
	var tmpQ = quat.create();
	var tmpV = vec3.create();
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [-dy, -dx, 0], rotAngle);
			//mat4.translate(renderer.scene.mvMat, renderer.scene.mvMat, [10*dx/canvas.width, -10*dy/canvas.height, 0]);

			//put tmpQ into the frame of renderer.view.angle, so we can rotate the view vector by it
			//  lastMouseRot = renderer.view.angle-1 * tmpQ * renderer.view.angle
			// now newViewAngle = renderer.view.angle * tmpQ = lastMouseRot * renderer.view.angle
			// therefore lastMouseRot is the global transform equivalent of the local transform of tmpQ
			quat.mul(lastMouseRot, renderer.view.angle, tmpQ);
			quat.conjugate(tmpQ, renderer.view.angle);
			quat.mul(lastMouseRot, lastMouseRot, tmpQ);

			vec3.copy(tmpV, renderer.view.pos);
			var posDist = vec3.length(tmpV);
			vec3.transformQuat(tmpV, tmpV, lastMouseRot);
			vec3.normalize(tmpV, tmpV);
			vec3.scale(tmpV, tmpV, posDist);
			vec3.copy(renderer.view.pos, tmpV);
			
			//RHS apply so it is relative to current view 
			//newViewAngle := renderer.view.angle * tmpQ
			quat.mul(renderer.view.angle, lastMouseRot, renderer.view.angle);
			quat.normalize(renderer.view.angle, renderer.view.angle);

		},
		zoom : function(dz) {
			renderer.view.fovY *= Math.exp(-.0003 * dz);
			renderer.view.fovY = Math.clamp(renderer.view.fovY, 1, 179);
			renderer.updateProjection();
		}
	});
	resize();
	update();
}


