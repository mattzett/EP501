%% Need some linear algebra tools in order to solve elliptic equations
addpath ../linear_algebra;


%% Define a 1D space and time grid in x for test problem
lx=64;
a=0;     %here a,b are the endpoints of the x-domain
b=1;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing


%% Define parameters of the parabolic equation, time variable
lambda=2;
tau=1/(2*pi/(2*dx))^2/lambda;    %diffusion time scale for the equation, based on smallest resolvable spatial mode
%dt=tau/5;              %time step

dtmargin=1/lambda/2*dx^2;
dt=25*dtmargin;
tmin=0;
tmax=1024*tau;          %go out to three times the diffusion time scale for the smallest possible mode
t=tmin:dt:tmax;
lt=numel(t);


%% FTCS implementation
f=zeros(lx,lt);
f(:,1)=sin(2*pi*x)+sin(8*pi*x);

%FTCS iterations
for n=1:lt-1
    f(1,n+1)=0;   %assume temperature goes to some small number on left boundary
    for i=2:lx-1     %interior grid points
        f(i,n+1)=f(i,n)+dt/dx^2*lambda*(f(i+1,n)-2*f(i,n)+f(i-1,n));
    end %for
    f(lx,n+1)=0;  %assume temperature goes to some small number on right boundary
end %for

figure(1);
subplot(131);
imagesc(t,x,f);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)')
title('FTCS')
set(gca,'FontSize',16);


% %% Creation of a Matlab movie
% figure(2);
% for n=1:lt
%     plot(x,f(:,n));
%     set(gca,'FontSize',16);
%     xlabel('x (m)');
%     ylabel('f(x)');
%     title(sprintf('f(x) @ t=%f',t(n)))
%     axis([0 1 -2 2]);
%     M(n)=getframe;
% end %for
% %movie(M);   %for whatever reason this doesn't store the axis labels and
% %title which makes it kind of worthless.  


%% Trapezoidal implementation, note matrix solutions are more efficiently handled thru tri-diagonal solver; Matlab built-in will detect automatically
f2=zeros(lx,lt);
A=sparse(lx,lx);   %allocate sparse array storage (this matrix is to be tridiag)
b=zeros(lx,1);


f2(:,1)=sin(2*pi*x)+sin(8*pi*x);
for n=2:lt
    A(1,1)=1;
    b(1)=0;
    for ix=2:lx-1
        %i-1 coeff
        A(ix,ix-1)=-lambda/2/dx^2;
        
        %i coeff
        A(ix,ix)=1/dt+lambda/dx^2;
        
        %i+1 coeff
        A(ix,ix+1)=-lambda/2/dx^2;
        
        b(ix)=f2(ix,n-1)/dt+(f2(ix+1,n-1)-2*f2(ix,n-1)+f2(ix-1,n-1))/dx^2*(lambda/2);
    end %for
    A(lx,lx)=1;
    b(lx)=0;
    
    fnow=A\b;
    f2(:,n)=fnow;
end %for



%% Compare two solutions on plot
figure(1);
subplot(132);
imagesc(t,x,f2);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('C-N solution');
set(gca,'FontSize',16);


%% Compute and plot the analytical solution (see https://github.com/gemini3d/GEMINI-docs/blob/master/test_descriptions/GEMINItests.pdf for a derivation)
[T,X]=meshgrid(t,x);
tempexact=exp(-4*pi^2*lambda*T).*sin(2*pi*X)+exp(-64*pi^2*lambda*T).*sin(8*pi*X);

figure(1);
subplot(133);
imagesc(t,x,tempexact);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('Exact');
set(gca,'FontSize',16);

%% Illustrate stability with FTCS


%% Adaptive time stepping?


%% Reset paths when we are done (for consistency, cleanliness)
rmpath ../linear_algebra;