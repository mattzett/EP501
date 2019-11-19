%% Define a 1D space and time grid in x,t for a test problem
lx=64;
N=lx*ly;
a=1;     %here a,b are the endpoints of the x-domain
b=a;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing
dt=0.1;              %time step


%% Lax method


%% Lax-Wendroff


%% Upwinding, Godunov's method


%% Illustration of conditional stability


