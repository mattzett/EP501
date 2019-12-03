%% Define a 1D space and time grid in x,t for a test problem
lx=64;
a=0;     %here a,b are the endpoints of the x-domain
b=1;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing
dt=0.015;              %time step
t=0:dt:10;
lt=numel(t);


%% Paramters of the PDE
v=1;            %velocity of wave propagation
CFL=v*dt/dx;    %CFL number should be <=1


%% Set the initial conditions
x0=1/2*(a+b);
sigx=1/15*(b-a);
finitial=exp(-(x-x0).^2/2/sigx^2);

figure(1);
plot(x,finitial);
xlabel('x');
ylabel('f(x,t_0)');
set(gca,'FontSize',24);


%% Analytical solution
% fexact=zeros(lx,lt);
% fexact(:,1)=finitial;
% 
% for n=1:lt-1
%     for i=1:lx
%       xloc=mod(x0+v*t(n),b);
%       fexact(i,n+1)=exp(-(x(i)-xloc).^2/2/sigx^2);
%     end %for
%     
%     plot(x,fexact(:,n));
%     xlabel('x');
%     ylabel('f(x,t)');
%     title('Exact')
%     pause(0.1);
% end %for


% %% FTCS method (will be unstable)
% figure(2);
% 
% fFTCS=zeros(lx,lt);
% fFTCS(:,1)=finitial;
% 
% %ghost cell values for implementing boundary conditions
% fleft=zeros(2,1);
% fright=zeros(2,1);
% 
% for n=1:lt-1
%     fleft=fFTCS(lx-1:lx,n);    %ghost cells here implement periodic boundary conditions
%     fright=fFTCS(1:2,n);
%     
%     fFTCS(1,n+1)=fFTCS(1,n)-dt/2/dx*v*(fFTCS(2,n)-fleft(2));
%     for i=2:lx-1     %interior grid points
%         fFTCS(i,n+1)=fFTCS(i,n)-dt/2/dx*v*(fFTCS(i+1,n)-fFTCS(i-1,n));
%     end %for
%     fFTCS(lx,n+1)=fFTCS(lx,n)-dt/2/dx*v*(fright(1)-fFTCS(lx-1,n));
%     
%     plot(x,fFTCS(:,n));
%     xlabel('x');
%     ylabel('f(x,t)');
%     title('FTCS');
%     set(gca,'FontSize',24);
%     pause(0.1);
% end %for
% 
% 
% %% Lax-Friedrich's method
% figure(3);
% 
% flax=zeros(lx,lt);
% flax(:,1)=finitial;
% 
% %ghost cell values for implementing boundary conditions
% fleft=zeros(2,1);
% fright=zeros(2,1);
% 
% for n=1:lt-1
%     fleft=flax(lx-1:lx,n);    %ghost cells here implement periodic boundary conditions
%     fright=flax(1:2,n);
%     
%     flax(1,n+1)=1/2*(fleft(2)+flax(2,n))-dt/2/dx*v*(flax(2,n)-fleft(2));
%     for i=2:lx-1     %interior grid points
%         flax(i,n+1)=1/2*(flax(i-1,n)+flax(i+1,n))-dt/2/dx*v*(flax(i+1,n)-flax(i-1,n));
%     end %for
%     flax(lx,n+1)=1/2*(flax(lx-1,n)+fright(1))-dt/2/dx*v*(fright(1)-flax(lx-1,n));
%     
%     plot(x,flax(:,n));
%     xlabel('x');
%     ylabel('f(x,t)');
%     title('L-F');
%     set(gca,'FontSize',24);
%     pause(0.1);
% end %for


%% Lax-Wendroff
figure(4);

flw=zeros(lx,lt);
flw(:,1)=finitial;

%ghost cell values for implementing boundary conditions
fleft=zeros(2,1);
fright=zeros(2,1);
fhalf=zeros(lx+1,1);

for n=1:lt-1
    fleft=flw(lx-1:lx,n);    %ghost cells here implement periodic boundary conditions
    fright=flw(1:2,n);
    
    %half step lax-f update for cell edges. note indexing for fhalf,
    %i->i-1/2, i+1->i+1/2
    fhalf(1)=1/2*(fleft(2)+flw(1,n))-dt/2/dx*v*(flw(1,n)-fleft(2));
    for i=2:lx     %interior grid points
        fhalf(i)=1/2*(flw(i-1,n)+flw(i,n))-dt/2/dx*v*(flw(i,n)-flw(i-1,n));
    end %for
    fhalf(lx+1)=1/2*(flw(lx,n)+fright(1))-dt/2/dx*v*(fright(1)-flw(lx,n));
    
    %full time step LW update
    flw(1,n+1)=flw(1,n)-dt/dx*v*(fhalf(2)-fhalf(1));
    for i=2:lx-1     %interior grid points    
      flw(i,n+1)=flw(i,n)-dt/dx*v*(fhalf(i+1)-fhalf(i));
    end %for
    flw(lx,n+1)=flw(lx,n)-dt/dx*v*(fhalf(lx+1)-fhalf(lx));
    
    plot(x,flw(:,n));
    xlabel('x');
    ylabel('f(x,t)');
    title('L-W');
    set(gca,'FontSize',24);
    pause(0.1);
end %for


%% Upwinding 


%% Problems with implicit approaches to waves


%% Godunov's method


