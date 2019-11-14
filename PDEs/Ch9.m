%% Need some linear algebra tools in order to solve elliptic equations
addpath ../linear_algebra;


%% Define a 2D grid in x,y for a test problem
lx=25;
ly=25;
N=lx*ly;
a=1;
b=a;     %use a square region for a test problem
x=linspace(0,a,lx);
y=linspace(0,b,ly);
dx=x(2)-x(1);    %constant grid spacing
dy=y(2)-y(1);    %ditto
[X,Y]=meshgrid(x,y);


%% Define Dirichlet boundary conditions for the test problem:
%https://github.com/gemini3d/GEMINI-docs/blob/master/test_descriptions/GEMINItests.pdf
f1=zeros(lx,1);
f2=sin(2*pi*x);
g1=zeros(1,ly);
g2=zeros(1,ly);


%% Setup of matrix for solving FDEs system size is NxN=lx*ly x lx*ly
M=zeros(N,N);
for j=1:ly
    for i=1:lx
        k=(j-1)*lx+i;
        if(j==1)      %min y
            M(k,k)=1;
            
            %RHS
            b(k)=f1(i);
        elseif(j==ly) %max y
            M(k,k)=1;
            
            %RHS
            b(k)=f2(i);            
        elseif(i==1)    %min x
            M(k,k)=1;
            
            %RHS
            b(k)=g1(j);           
        elseif(i==lx)   %max x
            M(k,k)=1;
            
            %RHS
            b(k)=g2(j);            
        else
          M(k,k-lx)=1/dy^2;
          M(k,k-1)=1/dx^2;
          M(k,k)=-2/dx^2-2/dy^2;
          M(k,k+1)=1/dx^2;
          M(k,k+lx)=1/dy^2;
          
          %RHS
          b(k)=0;
        end %if
    end %for
end %for



%% Execute of solution for our matrix system with self-coded solver from repo (direct)



%% Solution with the built-in Matlab matrix solver



%% Solution with Jacobi iterative solver from repo



%% Plot our solution



%% Compute and plot the analytical solution (see https://github.com/gemini3d/GEMINI-docs/blob/master/test_descriptions/GEMINItests.pdf for a derivation)



%% Reset paths when we are done (for consistency, cleanliness)
rmpath ../linear_algebra;

