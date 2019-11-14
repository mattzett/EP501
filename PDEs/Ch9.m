%% Need some linear algebra tools in order to solve elliptic equations
addpath ../linear_algebra;


%% Define a 2D grid in x,y for a test problem
lx=50;
ly=50;
a=1;
b=a;     %use a square region for a test problem
x=linspace(0,a,lx);
y=linspace(0,b,ly);
dx=x(2)-x(1);    %constant grid spacing
dy=y(2)=y(1);    %ditto
[X,Y]=meshgrid(x,y);


%% Define Dirichlet boundary conditions for the test problem



%% Execute of solution for our matrix system with self-coded solver from repo (direct)



%% Solution with the built-in Matlab matrix solver



%% Solution with Jacobi iterative solver from repo



%% Plot our solution



%% Compute and plot the analytical solution (see https://github.com/gemini3d/GEMINI-docs/blob/master/test_descriptions/GEMINItests.pdf for a derivation)



%% Reset paths when we are done (for consistency, cleanliness)
rmpath ../linear_algebra;

