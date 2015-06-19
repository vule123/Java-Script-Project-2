// constants
// max and min orbiting speed for planets
const maxSpeed = 3;
const minSpeed = 0.5;
// max radius of the sun
const maxSunRad = 12;
// complexity levels
const lowComplexity = 1;
const medComplexity = 2;
const highComplexity = 4;
// types of shading
const flatShading = 0;
const gouraudShading = 1;
const phongShading = 2;
// incremental adjustment of the camera's orientation
const deltaYaw = 1;
const deltaAngle = 3;
// incremental adjustment of the camera's horizontal view
const deltaView = 1;

// global variables
var planetsArray = [];
// different complexities for the four planets
var complexityLevelsArray = [lowComplexity, medComplexity, highComplexity, medComplexity];
var normalTypesArray = ['flat', 'smooth', 'smooth', 'smooth'];
var firstSurface = new Reflectance(vec3(0.25, 0.25, 0.25), vec3(0.4, 0.4, 0.4), vec3(0.774597, 0.774597, 0.774597), 76.8, flatShading); // icy white, faceted, diamond-like 
var secondSurface = new Reflectance(vec3(0.0215, 0.1745, 0.0215), vec3(0.07568, 0.61424, 0.07568), vec3(0.633, 0.727811, 0.633), 76.8, gouraudShading); // swampy, watery green 
var thirdSurface = new Reflectance(vec3(0.0, 0.05, 0.7), vec3(0.4, 0.5, 0.4), vec3(0.04, 0.7, 0.04), 10.0, phongShading);   // covered in clam smooth water 
var fourthSurface = new Reflectance(vec3(0.19125, 0.0735, 0.0225), vec3(0.7038, 0.27048, 0.0828), vec3(0.0, 0.0, 0.0), 1.0, phongShading); // covered in mud, brownish-orange

var reflectanceArray = [firstSurface, secondSurface, thirdSurface, fourthSurface];

var sun = new Star();
// normal camera
var camera = new Camera();
// attached camera
var attachedCamera = new Camera();
// attached projection
var attachedProjection = new Projection();
var isCameraAttached = false;
// normal projection
var normProjection = new Projection();

// Some supporting functions
// Scale a vector
function scaleVec(s, u, homogeneous) {
    var result = [];
    for (var i = 0; i < u.length - 1; ++i) 
        result.push(s * u[i]);

    if (homogeneous) 
        result.push(u[u.length - 1]);
    else
        result.push(s * u[u.length - 1]);

    return result;
}
// Multiplication of a matrix and a vector
function multMV(mat, vec) {
    if (mat.length == 0) {
        throw ("Matrix Vector Multiply: dimension wrong!");
    }
    for (var i = 0; i < mat.length; i++) {
        if (mat[i].length != vec.length) {
            throw ("Matrix Vector Multiply: dimension wrong!");
        }
    }

    var result = [];
    for (var i = 0; i < mat.length; i++) {
        var sum = 0.0;
        for (var j = 0; j < vec.length; j++) {
            sum += mat[i][j] * vec[j];
        }
        result.push(sum);
    }

    return result.splice(0, mat.length);
}
// return the 3x3 submatrix on the up left of a 4x4 matrix
function subMat(mat) {
    var result = mat3();
    for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) result[i][j] = mat[i][j];
    }
    return result;
}
// inverse of a 3x3 matrix
function mat3Invert(mat) {
    var result = mat3();
    result[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    result[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
    result[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    result[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    result[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
    result[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
    result[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
    result[2][1] = mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1];
    result[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    var det = mat[0][0] * result[0][0] + mat[0][1] * result[1][0] + mat[0][2] * result[2][0];
    if (det == 0) return mat3(vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0));
    for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) result[i][j] /= det;
    }
    return result;
}
// perspective projection
function perspective(near, far, top, bottom, left, right) {
    var r00 = 2 * near / (right - left);
    var r02 = (right + left) / (right - left);
    var r11 = 2 * near / (top - bottom);
    var r12 = (top + bottom) / (top - bottom);
    var r22 = -(far + near) / (far - near);
    var r23 = -2 * far * near / (far - near);

    var result = mat4();
    result[0][0] = r00;
    result[0][2] = r02;
    result[1][1] = r11;
    result[1][2] = r12;
    result[2][2] = r22;
    result[2][3] = r23;
    result[3][2] = -1;
    result[3][3] = 0;

    return result;
}

// Definitions of all the objects

function Planet() {
    this.pointsArray = [];  	// array of vertices forming a planet
    this.normalArray = [];		// normals at certain points in the pointsArray
    this.reflectance;			// surface of the planet
    this.satellitesArray = [];	// satellites orbiting the planet 

    this.centerPosArray = [];	// store the center of primitives for flat shading

    // calculate vertices for a unit sphere in the object coordinates
    this.tetrahedron = function (a, b, c, d, n, normalType) {
        this.divideTriangle(a, b, c, n, normalType);
        this.divideTriangle(d, c, b, n, normalType);
        this.divideTriangle(a, d, b, n, normalType);
        this.divideTriangle(a, c, d, n, normalType);
    }

    this.divideTriangle = function (a, b, c, count, normalType) {
        if (count <= 0) {
            this.triangle(a, b, c, normalType);
        }
        else {
            var ab = normalize(mix(a, b, 0.5), true);
            var ac = normalize(mix(a, c, 0.5), true);
            var bc = normalize(mix(b, c, 0.5), true);
            this.divideTriangle(a, ab, ac, count - 1, normalType);
            this.divideTriangle(ab, b, bc, count - 1, normalType);
            this.divideTriangle(bc, c, ac, count - 1, normalType);
            this.divideTriangle(ab, bc, ac, count - 1, normalType);
        }
    }

    this.triangle = function (a, b, c, normalType) {
        this.pointsArray.push(a);
        this.pointsArray.push(b);
        this.pointsArray.push(c);

        // generate different types of normal
        if (normalType == 'smooth') {
            // use the current normal of the points for 'smooth' normal type
            this.normalArray.push(a);
            this.normalArray.push(b);
            this.normalArray.push(c);
        }
        else if (normalType == 'flat') {
            // calculate the normal 
            var ab = subtract(b, a);
            var ac = subtract(c, a);
            var normal = normalize(vec4(cross(ac, ab), 0), false);
            // assign them to the three vertices
            this.normalArray.push(normal);
            this.normalArray.push(normal);
            this.normalArray.push(normal);
        }

    }

    this.transformation = mat4();			// the object-world transformation of this object


    this.init = function (rad, center, complexity, normalType, reflectance) {
        // generate a unit sphere in the object coordinate
        var va = vec4(0.0, 0.0, -1.0, 1);
        var vb = vec4(0.0, 0.942809, 0.333333, 1);
        var vc = vec4(-0.816497, -0.471405, 0.333333, 1);
        var vd = vec4(0.816497, -0.471405, 0.333333, 1);

        this.tetrahedron(va, vb, vc, vd, complexity, normalType);

        if (reflectance.shade() == flatShading) {
            for (var i = 0; i < this.pointsArray.length; i += 3) {
                var centerPos = vec4();
                centerPos = add(centerPos, this.pointsArray[i]);
                centerPos = add(centerPos, this.pointsArray[i + 1]);
                centerPos = add(centerPos, this.pointsArray[i + 2]);
                centerPos = scaleVec(1 / 3, centerPos, false);
                centerPos[3] = 1;
                this.centerPosArray.push(centerPos);
                this.centerPosArray.push(centerPos);
                this.centerPosArray.push(centerPos);
            }
        }

        // generate the planet by scaling the unit sphere in the object coordinate
        for (var i = 0; i < this.pointsArray.length; i++) {
            this.pointsArray[i] = scaleVec(rad, this.pointsArray[i], true);
        }
        
        // transform the planet into the world coordinate
        this.translation(center[0], center[1], center[2]);

        // set the surface of the planet
        this.reflectance = reflectance;
    }
    this.attachSatellite = function (orbit) {
        this.satellitesArray.push(orbit);
    }
    this.translation = function (x, y, z) {
        this.transformation = mult(translate(x, y, z), this.transformation);

        // translate the entire system by moving all the planets and their satellites
        for (var i = 0; i < this.satellitesArray.length; i++) {
            this.satellitesArray[i].translation(x, y, z);
        }
    }
    this.rotation = function (theta, axis) {
        this.transformation = mult(rotate(theta, axis), this.transformation);

        for (var i = 0; i < this.satellitesArray.length; i++) {
            this.satellitesArray[i].rotation(theta, axis);
        }
    }
    this.scale = function (x, y, z) {
        this.transformation = mult(scale(x, y, z), this.transformation);
    }
    this.bind = function (gl) {
        // get the pipeline
        var program = gl.getParameter(gl.CURRENT_PROGRAM);

        // calculate model-view transformation by multiplying object-world and world-eye transformation

        var modelViewMatrix;
        if (isCameraAttached) modelViewMatrix = mult(attachedCamera.worldToEye(), this.transformation);
        else modelViewMatrix = mult(camera.worldToEye(), this.transformation);
        // calculate normal Matrix as the transpose of inverse of the up-left 3x3 submatrix of model-view transformation
        var normalMatrix = mat3Invert(subMat(modelViewMatrix));
        normalMatrix = transpose(normalMatrix);

        var modelLoc = gl.getUniformLocation(program, "modelViewMatrix");
        gl.uniformMatrix4fv(modelLoc, false, new Float32Array(flatten(modelViewMatrix)));

        var normMatLoc = gl.getUniformLocation(program, "normalMatrix");
        gl.uniformMatrix3fv(normMatLoc, false, new Float32Array(flatten(normalMatrix)));

        // bind vertex position
        var positionLoc = gl.getAttribLocation(program, "vPosition");
        var positionBuffer = gl.createBuffer()
        gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(flatten(this.pointsArray)), gl.DYNAMIC_DRAW);

        gl.enableVertexAttribArray(positionLoc);
        gl.vertexAttribPointer(positionLoc, 4, gl.FLOAT, false, 0, 0);

        // bind vertex normal
        var normalLoc = gl.getAttribLocation(program, "vNormal");
        var normalBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, normalBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(flatten(this.normalArray)), gl.DYNAMIC_DRAW);

        gl.enableVertexAttribArray(normalLoc);
        gl.vertexAttribPointer(normalLoc, 4, gl.FLOAT, false, 0, 0);

        // for flat shading, bind primitive center positions for all its vertexs
        if (this.reflectance.shade() == flatShading) {
            var centerPosLoc = gl.getAttribLocation(program, "vFlatPosition");
            var centerPosBuffer = gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER, centerPosBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(flatten(this.centerPosArray)), gl.DYNAMIC_DRAW);

            gl.enableVertexAttribArray(centerPosLoc);
            gl.vertexAttribPointer(centerPosLoc, 4, gl.FLOAT, false, 0, 0);
        }
        else {
            var centerPosLoc = gl.getAttribLocation(program, "vFlatPosition");
            gl.disableVertexAttribArray(centerPosLoc);
        }

        // bind the reflectance of the planet
        var ambientLoc = gl.getUniformLocation(program, "ambientProduct");
        gl.uniform3fv(ambientLoc, this.reflectance.ambientProduct(sun));

        var diffuseLoc = gl.getUniformLocation(program, "diffuseProduct");
        gl.uniform3fv(diffuseLoc, this.reflectance.diffuseProduct(sun));

        var specularLoc = gl.getUniformLocation(program, "specularProduct");
        gl.uniform3fv(specularLoc, this.reflectance.specularProduct(sun));

        var shinLoc = gl.getUniformLocation(program, "shininess");
        gl.uniform1f(shinLoc, this.reflectance.shine());

        // bind the shading type of the planet: flat, Gouraud or Phong
        var shadeLoc = gl.getUniformLocation(program, "shadingType");
        gl.uniform1i(shadeLoc, this.reflectance.shade());
    };
    this.draw = function (gl) {
        gl.drawArrays(gl.TRIANGLES, 0, this.pointsArray.length);

        // draw the whole galaxy of the planet by drawing all of its satellites
        for (var i = 0; i < this.satellitesArray.length; i++) {
            this.satellitesArray[i].rotate();
            this.satellitesArray[i].bind(gl);
            this.satellitesArray[i].draw(gl);
        }
    };
}

function Orbit(radius, center, planet) {
    this.radius = radius;
    this.center = center;
    this.planet = planet;
    this.cameraAttached = false;
    this.camera;
    // generate rotation speed of the planets randomly
    this.rotateSpeed = Math.random() * (maxSpeed - minSpeed) + minSpeed;
    
    this.attach = function (attached) {
        this.camera = attached;
        this.cameraAttached = true;
    }
    this.rotate = function () {
        // rotate the planets around the center of the orbit in the world coordinate
        this.planet.translation(-this.center[0], -this.center[1], -this.center[2]);
        this.planet.rotation(this.rotateSpeed, [0, 1, 0]);
        this.planet.translation(this.center[0], this.center[1], this.center[2]);

        if (this.cameraAttached) {
            var delta = subtract(this.center.slice(0, 3), this.camera.position);
            delta = this.camera.deltaModify(vec4(delta, 0));
            this.camera.translation(delta[0], delta[1], delta[2]);
            this.camera.rotation(this.rotateSpeed, [0, 1, 0]);
            this.camera.translation(-delta[0], -delta[1], -delta[2]);
        }
    }
    this.rotation = function (theta, axis) {
        // rotate the entire orbit
        this.center = multMV(rotate(theta, axis), this.center);
        this.planet.rotation(theta, axis);

        if (this.cameraAttached) {
            this.camera.rotation(theta, axis);
        }
    }
    this.translation = function (x, y, z) {
        // translate the entire orbit
        this.center = multMV(translate(x, y, z), this.center);
        this.planet.translation(x, y, z);

        if (this.cameraAttached) {
            this.camera.translation(x, y, z);
        }
    }
    this.bind = function (gl) {
        this.planet.bind(gl);
    }
    this.draw = function (gl) {
        this.planet.draw(gl);
    }
}

function Reflectance(ambient, diffuse, specular, shininess, shadeType) {
    this.ambient = ambient;
    this.diffuse = diffuse;
    this.specular = specular;
    this.shininess = shininess;
    this.shadeType = shadeType;
    
    this.ambientProduct = function (light) {
        return mult(this.ambient, light.ambient());
    }
    this.diffuseProduct = function (light) {
        return mult(this.diffuse, light.diffuse());
    }
    this.specularProduct = function (light) {
        return mult(this.specular, light.specular());
    }
    this.shine = function () {
        return this.shininess;
    }
    this.shade = function () {
        return this.shadeType;
    }
}

// A star is at the center of the system (EX: the sun)
function Star() {
    this.planet = new Planet();
    this.position = vec4();
    this.lightAmbient = vec3();
    this.lightDiffuse = vec3();
    this.lightSpecular = vec3();
    
    this.init = function (rad, center, complexity, normalType, reflectance) {
        this.planet.init(rad, center, complexity, normalType, reflectance);

        // negate the normal of the star because it contains the light source
        for (var i = 0; i < this.planet.normalArray.length; i++) {
            this.planet.normalArray[i] = negate(this.planet.normalArray[i]);
        }
        // set up the light source
        this.position = center;
        this.lightAmbient = vec3(0.2, 0.2, 0.2);
        this.lightDiffuse = vec3(1.0, 1.0, 1.0);
        this.lightSpecular = vec3(1.0, 1.0, 1.0);
    }
    this.bind = function (gl) {
        var program = gl.getParameter(gl.CURRENT_PROGRAM);
        var lightLoc = gl.getUniformLocation(program, "lightPosition");

        // transform the light position from world coordinate to eye coordinate
        var lightPosition;
        if (isCameraAttached) lightPosition = multMV(attachedCamera.worldToEye(), this.position);
        else lightPosition = multMV(camera.worldToEye(), this.position);
        gl.uniform4fv(lightLoc, lightPosition);

        this.planet.bind(gl);
    }
    this.draw = function (gl) {
        this.planet.draw(gl);
    }
    this.ambient = function () {
        return this.lightAmbient;
    }
    this.diffuse = function () {
        return this.lightDiffuse;
    }
    this.specular = function () {
        return this.lightSpecular;
    }
}

function Camera() {
    this.initial = [vec3(0, 0, 0), vec3(0, 0, 0)];    // inital position and orientation
    this.position = vec3(0, 0, 0);   				  // current position
    this.orientation = vec4(0, 0, 1, 0);              // current orientation			
    this.top = vec4(0, 1, 0, 0);					  // current top orientation
    this.transformation = mat4();					  // world-eye transformation
    this.modification = mat4();						  // modification matrix is for converting delta vector between
                                                         // rotation center of the camera and the position of the camera
                                                         // from world coordinate to eye coordinate 
    
    this.init = function (x, y, z, theta) {
        this.initial = [vec3(x, y, z), theta];
        this.position = vec3(x, y, z);
        this.transformation = translate(-x, -y, -z);
        this.transformation = mult(rotate(-theta[0], [1, 0, 0]), this.transformation);
        this.transformation = mult(rotate(-theta[1], [0, 1, 0]), this.transformation);
        this.transformation = mult(rotate(-theta[2], [0, 0, 1]), this.transformation);

        // calculate the initial top and head orientation
        var inverseTransform = mat4();
        inverseTransform = mult(rotate(theta[0], [1, 0, 0]), inverseTransform);
        inverseTransform = mult(rotate(theta[1], [0, 1, 0]), inverseTransform);
        inverseTransform = mult(rotate(theta[2], [0, 0, 1]), inverseTransform);
        this.orientation = multMV(inverseTransform, vec4(0, 0, 1, 0));
        this.top = multMV(inverseTransform, vec4(0, 1, 0, 0));
    }
    this.translation = function (x, y, z) {
        this.transformation = mult(translate(-x, -y, -z), this.transformation);
        this.position = add(this.position, vec3(x, y, z));
    }
    this.rotation = function (theta, axis) {
        this.transformation = mult(rotate(-theta, axis), this.transformation);
        this.orientation = mult(rotate(theta, axis), this.orientation);
    }
    this.worldToEye = function () {
        return this.transformation;
    }
    // control the camera in the eye coordinate
    this.up = function () {
        var direction = scaleVec(0.25, this.top, false);
        this.transformation = mult(translate(-direction[0], -direction[1], -direction[2]), this.transformation);
        this.position = add(this.position, direction.slice(0, 3));
    }
    this.down = function () {
        var direction = scaleVec(-0.25, this.top, false);
        this.transformation = mult(translate(-direction[0], -direction[1], -direction[2]), this.transformation);
        this.position = add(this.position, direction.slice(0, 3));
    }
    this.forward = function () {
        var direction = scaleVec(0.25, this.orientation, false);
        this.transformation = mult(translate(direction[0], direction[1], direction[2]), this.transformation);
        this.position = add(this.position, negate(direction).slice(0, 3));
    }
    this.backward = function () {
        var direction = scaleVec(-0.25, this.orientation, false);
        this.transformation = mult(translate(direction[0], direction[1], direction[2]), this.transformation);
        this.position = add(this.position, negate(direction).slice(0, 3));
    }
    this.left = function () {
        var direction = cross(this.orientation, this.top);
        direction = normalize(direction, false);
        direction = scaleVec(-0.25, direction, false);
        this.transformation = mult(translate(direction[0], direction[1], direction[2]), this.transformation);
        this.position = add(this.position, negate(direction).slice(0, 3));
    }
    this.right = function () {
        var direction = cross(this.orientation, this.top);
        direction = normalize(direction, false);
        direction = scaleVec(0.25, direction, false);
        this.transformation = mult(translate(direction[0], direction[1], direction[2]), this.transformation);
        this.position = add(this.position, negate(direction).slice(0, 3));
    }
    this.lYaw = function () {
        // rotate the camera in the eye coordinate
        this.transformation = mult(rotate(-deltaYaw, this.top.slice(0, 3)), this.transformation);
        this.modification = mult(rotate(-deltaYaw, this.top.slice(0, 3)), this.modification);
        this.orientation = multMV(rotate(deltaYaw, this.top.slice(0, 3)), this.orientation);
    }
    this.rYaw = function () {
        this.transformation = mult(rotate(deltaYaw, this.top.slice(0, 3)), this.transformation);
        this.modification = mult(rotate(deltaYaw, this.top.slice(0, 3)), this.modification);
        this.orientation = multMV(rotate(-deltaYaw, this.top.slice(0, 3)), this.orientation);
    }
    this.initialStatus = function () {
        return this.initial;
    }
    this.deltaModify = function (delta) {
        return multMV(this.modification, delta).slice(0, 3);
    }
}

function Projection() {
    // the parameters of view box
    this.near = -1;
    this.far = 1;
    this.top = 1;
    this.bottom = -1;
    this.left = -1;
    this.right = 1;
    this.type = "orthographic";
    this.init = function (near, far, top, bottom, left, right, type) {
        this.near = near;
        this.far = far;
        this.top = top;
        this.bottom = bottom;
        this.left = left;
        this.right = right;
        this.type = type;
    }
    this.bind = function (gl) {
        var program = gl.getParameter(gl.CURRENT_PROGRAM);

        // calculate the projection matrix 
        var projection = mat4();
        if (this.type == "perspective") {
            projection = perspective(this.near, this.far, this.top, this.bottom, this.left, this.right);

            var projLoc = gl.getUniformLocation(program, "perspectiveProjection");
        }
        else if (this.type == "orthographic") {
            projection = ortho(this.left, this.right, this.bottom, this.top, this.near, this.far);

            var projLoc = gl.getUniformLocation(program, "orthographicProjection");
        }

        gl.uniformMatrix4fv(projLoc, false, new Float32Array(flatten(projection)));
    }
    this.change = function (gl, type) {
        var program = gl.getParameter(gl.CURRENT_PROGRAM);

        var boolLoc = gl.getUniformLocation(program, "isPerspective");

        if (type == "perspective") {
            gl.uniform1i(boolLoc, true);
        }
        else if (type = "orthographic") {
            gl.uniform1i(boolLoc, false);
        }
    }
}

function init() {
    var canvas = document.getElementById("glcanvas");
    var gl = initWebGL(canvas);

    if (gl) {
        // setup the render pipeline consisting of programmable shaders
        var vertexShader = getShader(gl, "vertex-shader");
        var fragmentShader = getShader(gl, "fragment-shader");

        if (!vertexShader || !fragmentShader) {
            return;
        }

        var renderPipeline = getProgram(gl, [vertexShader, fragmentShader]);

        if (!renderPipeline) {
            return;
        }

        gl.useProgram(renderPipeline);

        // initiate sun
        // randomly generate radius of sun and calculate the color of sun according to its radius
        var sunRad = Math.random() * 6 + 4;
        // set the center of sun at (0, 0, -10, 1) in the world coordinate
        var sunCenter = vec4(0, 0, -10, 1);
        var sunReflect = sunRad / maxSunRad;
        var sunReflectance = new Reflectance(vec3(sunReflect, 0, 1 - sunReflect), vec3(sunReflect, 0, 1 - sunReflect), vec3(sunReflect, 0, 1 - sunReflect), 38.4, phongShading);
        sun.init(sunRad, sunCenter, highComplexity, 'smooth', sunReflectance);

        var systemRad = 0;
        // generate planets and orbits
        var rad, orbitRad, pos;
        for (var i = 0; i < complexityLevelsArray.length; i++) {
            planetsArray.push(new Planet());
            // randomly generate radius and position of the four planets and radius of their orbits
            rad = Math.random() * 2 + 1;
            orbitRad = Math.random() * 4 + 6 + 10 * (i + 1);

            // increase the radius of the orbit to allocate space for the moon
            if (i == complexityLevelsArray.length - 1) orbitRad += 2;

            // translate the center of the orbit from the origin to the center of the sun in the world coordinate
            pos = scaleVec(orbitRad, normalize(vec4(Math.random(), 0, Math.random(), 1), true), true);
            var fPos = add(pos, sunCenter);
            fPos[3] = 1;
            planetsArray[i].init(rad, fPos, complexityLevelsArray[i], normalTypesArray[i], reflectanceArray[i]);
            // attach camera to the second planet
            var tempOrbit = new Orbit(orbitRad, sunCenter, planetsArray[i]);
            if (i == 1) {
                fPos = add(fPos, vec4(0, rad, 0, 1));
                attachedCamera.init(fPos[0], fPos[1], fPos[2], vec3(0, 0, 0));
                tempOrbit.attach(attachedCamera);
            }
            // put the four planets around the sun
            sun.planet.attachSatellite(tempOrbit);
        }

        // generate moon and assign it to the last planet as its satellite
        var moon = new Planet();
        var moonRad = 0.5;
        var moonOrbitRad = rad + 1.5;
        var moonPos = scaleVec(orbitRad + moonOrbitRad, normalize(pos, true), true);
        var moonReflectance = new Reflectance(vec3(0.25, 0.25, 0.25), vec3(0.4, 0.4, 0.4), vec3(0.774597, 0.774597, 0.774597), 76.8, flatShading);
        pos = scaleVec(orbitRad, pos, true);
        var fPos = add(moonPos, sunCenter);
        fPos[3] = 1;
        moon.init(moonRad, fPos, medComplexity, 'flat', moonReflectance);
        fPos = add(pos, sunCenter);
        fPos[3] = 1;
        var moonOrbit = new Orbit(moonOrbitRad, fPos, moon);
        planetsArray[i - 1].attachSatellite(moonOrbit);

        systemRad = orbitRad + moonOrbitRad + moonRad;

        // initial model-view transformation to place the camera
        camera.init(0, systemRad, 2 * systemRad, vec3(-30, 0, 0));

        // calculate perspective projection
        // bind the perspective projection matrix
        normProjection.init(systemRad + 20, 3 * systemRad + 20, systemRad, -systemRad, -systemRad, systemRad, "perspective");
        normProjection.bind(gl);
        // bind the orthographic projection matrix
        normProjection.init(systemRad + 20, 3 * systemRad + 20, systemRad, -systemRad, -systemRad, systemRad, "orthographic");
        normProjection.bind(gl);
        // change the projection type
        normProjection.change(gl, "perspective");
        // initiate the attached camera
        attachedProjection.init(-sunRad, systemRad, sunRad, -sunRad, -sunRad, sunRad, "orthographic");


        gl.enable(gl.DEPTH_TEST);
        gl.clearDepth(1.0, 1.0, 1.0, 1.0);
        gl.clearColor(0.0, 0.0, 0.0, 1.0);
        gl.clear(gl.DEPTH_BUFFER_BIT);
        render();
    }
}

function initWebGL(canvas) {
    var gl = null;

    // set up the canvas in webGL, and check if the browser supports it
    try {
        gl = canvas.getContext("webgl") || canvas.getContext("experimental-webgl");
    }
    catch (e) { }

    if (!gl) {
        alert("Unable to initialize WebGL. Your browser may not support it.");
        gl = null;
    }
    return gl;
}

function getShader(gl, id) {
    var shaderScript = document.getElementById(id);

    if (!shaderScript) {
        console.log("Unable to find shader.");
        return null;
    }

    // get the GLSL source code of shader
    var Source = "";
    currentChild = shaderScript.firstChild;
    while (currentChild) {
        if (currentChild.nodeType == currentChild.TEXT_NODE) {
            Source += currentChild.textContent;
        }
        currentChild = currentChild.nextSibling;
    }

    // create shader of the given type
    var Shader;
    if (shaderScript.type == "x-shader/x-vertex") {
        Shader = gl.createShader(gl.VERTEX_SHADER);
    }
    else if (shaderScript.type == "x-shader/x-fragment") {
        Shader = gl.createShader(gl.FRAGMENT_SHADER);
    }
    else {
        console.log("Incorrect type of shader.");
        return null;
    }

    // check if the compilation of the shader is successful
    gl.shaderSource(Shader, Source);
    gl.compileShader(Shader);

    if (!gl.getShaderParameter(Shader, gl.COMPILE_STATUS)) {
        console.log("Shader compilation error: " + gl.getShaderInfoLog(Shader));
        return null;
    }

    return Shader;
}

function getProgram(gl, arrayShader) {
    var Program = gl.createProgram();

    // attach shaders to the program
    for (var i = 0; i < arrayShader.length; i++) {
        gl.attachShader(Program, arrayShader[i]);
    }

    gl.linkProgram(Program);

    if (!gl.getProgramParameter(Program, gl.LINK_STATUS)) {
        console.log("Program link error: " + gl.getProgramInfoLog(Program));
        return null;
    }

    return Program;
}

function render() {
    setTimeout(function () {
        window.requestAnimationFrame(render);
        var gl = document.getElementById("glcanvas").getContext("experimental-webgl");
        gl.clear(gl.COLOR_BUFFER_BIT);
        sun.bind(gl);
        sun.draw(gl);

    }, 16)
}

// Use keyboard to manipulate movement of the camera
function keyPressed(event) {
    var evtObject = window.event ? event : e; // event for IE, e for FireFox
    var keyUnicode = evtObject.charCode ? evtObject.charCode : evtObject.keyCode;
    var Key = String.fromCharCode(keyUnicode);

    // get the canvas of WebGL to implement events
    var gl = document.getElementById("glcanvas").getContext("experimental-webgl");
    if (!gl) {
        return null;
    }

    if (Key == '&') { // Up
        if (!isCameraAttached) {
            camera.up();
        }
    }
    else if (Key == '(') { // Down
        if (!isCameraAttached) {
            camera.down();
        }
    }
    else if (Key == 'I') { // Forward
        if (!isCameraAttached) {
            camera.forward();
        }
    }
    else if (Key == 'M') { // Backward
        if (!isCameraAttached) {
            camera.backward();
        }
    }
    else if (Key == 'J') { // Left
        if (!isCameraAttached) {
            camera.left();
        }
    }
    else if (Key == 'K') { // Right
        if (!isCameraAttached) {
            camera.right();
        }
    }
    else if (Key == '\'') { // Right heading
        if (!isCameraAttached) camera.rYaw();
        else attachedCamera.rYaw();
    }
    else if (Key == '%') { // Left heading
        if (!isCameraAttached) camera.lYaw();
        else attachedCamera.lYaw();
    }
    else if (Key == 'B') { // reset to the inital position
        if (!isCameraAttached) {
            var initial = camera.initialStatus();
            camera.init(initial[0][0], initial[0][1], initial[0][2], initial[1]);
        }
    }
    else if (Key == 'A') { // attach/detach the camera
        if (isCameraAttached) {
            isCameraAttached = false;
            normProjection.bind(gl);
            normProjection.change(gl, "perspective");
        }
        else {
            isCameraAttached = true;
            attachedProjection.bind(gl);
            attachedProjection.change(gl, "orthographic");
        }
    }
    else if (Key == 'N') { // narrow the horizontal view 
        if (normProjection.left + deltaView / 2 < normProjection.right - deltaView / 2 && normProjection.bottom + deltaView / 2 < normProjection.top - deltaView / 2) {
            normProjection.init(normProjection.near, normProjection.far, normProjection.top - deltaView / 2, normProjection.bottom + deltaView / 2, normProjection.left + deltaView / 2, normProjection.right - deltaView / 2, "perspective");
        }
        normProjection.bind(gl);
    }
    else if (Key == 'W') { // widen the horizontal view 
        normProjection.init(normProjection.near, normProjection.far, normProjection.top + deltaView / 2, normProjection.bottom - deltaView / 2, normProjection.left - deltaView / 2, normProjection.right + deltaView / 2, "perspective");
        normProjection.bind(gl);
    }
}