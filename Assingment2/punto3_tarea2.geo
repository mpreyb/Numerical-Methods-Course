// Gmsh project created on Thu Oct 07 23:44:06 2021
r_in = 1.0;   // Internal radius
r_ext = 1.2;  // External radius


// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {r_in, 0, 0, 1.0};
Point(3) = {-r_in, 0, 0, 1.0};
Point(4) = {r_ext, 0, 0, 1.0};
Point(5) = {-r_ext, 0, 0, 1.0};

// Circulos
Circle(1) = {5, 1, 4};
Circle(2) = {5, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {3, 1, 2};
Circle(5) = {2, 1, 3};

//Superficies
Curve Loop(1) = {5, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {3, 1};
Plane Surface(2) = {2};

// Grupo físico
Physical Curve(6) = {5, 4}; //interior
Physical Curve(7) = {3, 1}; //exterior
Physical Surface(8) = {1}; //adentro
Physical Surface(9) = {2}; //afuera
