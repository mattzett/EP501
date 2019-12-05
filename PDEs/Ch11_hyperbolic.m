%% Define a 1D space and time grid in x,t for a test problem
lx=64;
a=0;     %here a,b are the endpoints of the x-domain
b=1;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing

v=1;            %velocity of wave propagation
targetCFL=0.9;
dt=targetCFL*dx/v;
N=75;           %number of time steps to take
t=0:dt:N*dt;
lt=numel(t);


%% Initial conditions for our test problem
x0=1/2*(a+b);
sigx=1/15*(b-a);
%finitial=exp(-(x-x0).^2/2/sigx^2);
finitial=exp(-(x-x0).^20/2/sigx^20);


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


%% Lax-Friedrichs, Lax-Wendroff, Upwind (Godunov), and BTCS comparison
figure(3);

flax=zeros(lx,lt);
flax(:,1)=finitial;

flw=zeros(lx,lt);
flw(:,1)=finitial;

fgod=zeros(lx,lt);
fgod(:,1)=finitial;

fBTCS=zeros(lx,lt);
fBTCS(:,1)=finitial;
for n=1:lt-1
    flax(:,n+1)=LaxFried(dt,dx,v,flax(:,n));
    flw(:,n+1)=LaxWen(dt,dx,v,flw(:,n));
    fgod(:,n+1)=Godunov(dt,dx,v,fgod(:,n));
    fBTCS(:,n+1)=BTCS_hyp(dt,dx,v,fBTCS(:,n));
        
    plot(x,flax(:,n+1),x,flw(:,n+1),x,fgod(:,n+1),x,fBTCS(:,n+1));
    legend('Lax-F','Law-W','Upwind','BTCS');
    xlabel('x');
    ylabel('f(x,t)');
    title(sprintf('Solver comparison, t=%5.3f',t(n)));
    set(gca,'FontSize',24);
    
    if (n==1)
        pause;
    else
        pause(0.1);
    end %if
end %for
