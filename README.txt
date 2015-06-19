Name: Vu Le
UID: 004497690

Assignment 2:

1. Display a WebGL capable HTML canvas of size 960x540. The canvas has the z-buffer enabled and cleared to a black background.

2. Implement the function to generate a unit sphere, based on the function in the textbook, with a parameter to define the number of vertices that form the sphere.

3. Use the parameter ‘normalType’ to control the type of normal vector (flat or smooth). Normal vectors are generated in association with corresponding vertices.

4. Create a solar system with one sun and four orbiting planets. The center of the solar system is at (0, 0, -10) in the world coordinate. The radius of sun is random, ranging from 4 to 12. The radiuses of the four planets are random, ranging from 1 to 3. The orbiting speed of the planets is also random.

5. The sun is a light source. The sun’s size relative to the planets determines its color. The larger it is, the warmer its color is. The appearance of each planet is based on the descriptions. 

6. The flat and Gouraud shading are calculated in the vertex shader. The Phong shading is calculated in the fragment shader. The normal matrices, normal vectors, the positions of centers of primitives in the object coordinate, and the light’s position in the eye coordinate are calculated in the CPU.

7. Re-use the keyboard based navigation system from assignment 1 to allow moving around the solar system. Key ‘B’ is to reset the view of the system and its initial position.


EXTRA CREDIT

1. A moon is added to the last planet.

2. Pressing the key ‘A’ is to attach/detach the camera to the second planet. 
