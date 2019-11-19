%% Need some linear algebra tools in order to solve elliptic equations
addpath ../linear_algebra;


%% Define a 1D space and time grid in x,t for a test problem
lx=25;
ly=25;
N=lx*ly;
a=1;     %here a,b are the endpoints of the x-domain
b=a;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing
dt=0.1;              %time step


%% Parameters of the parabolic equation
lambda=2;


%% FTCS implementation


%% Trapezoidal implementation


%% Compare two solutions


%% Compute and plot the analytical solution (see https://github.com/gemini3d/GEMINI-docs/blob/master/test_descriptions/GEMINItests.pdf for a derivation)


%% Reset paths when we are done (for consistency, cleanliness)
rmpath ../linear_algebra;