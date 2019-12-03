clear,clc,close all;


%% Gridding in time
N=25;
tmin=0;
tmax=10;
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);


%% Test problem
y0=1;
alpha=2;
ybar=y0*exp(-alpha*t);
figure(1);
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);


%% Forward Euler solution
yfwd=zeros(1,N);
yfwd(1)=y0;
for n=2:N
    yfwd(n)=yfwd(n-1)*(1-alpha*dt);
end %for
hold on;
plot(t,yfwd,'--');


%% Backward Euler solution
ybwd=zeros(1,N);
ybwd(1)=y0;
for n=2:N
    ybwd(n)=ybwd(n-1)/(1+alpha*dt);
end %for
plot(t,ybwd,'-.');


%% Illustration of stability
alpha2=10;
dtmargin=2/alpha2;
dtstable=1/10*dtmargin;
dtunstable=1.01*dtmargin;

dt2=dtunstable;
t2=0:dt2:5;
N=numel(t);
y2bar=y0*exp(-alpha2*t2);

y2fwd=zeros(1,N);
y2fwd(1)=y0;
for n=2:N
    y2fwd(n)=y2fwd(n-1)*(1-alpha2*dt2);
end %for

figure;
plot(t2,y2bar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
hold on;
plot(t2,y2fwd,'--');


%% Illustration of convergence
%run the code above with finer mesh, etc.


%% Second order method; RK2
yRK2=zeros(1,N);
yRK2(1)=y0;
for n=2:N
    yhalf=yRK2(n-1)+dt/2*(-alpha*yRK2(n-1));
    yRK2(n)=yRK2(n-1)+dt*(-alpha*yhalf);
end %for
figure(1);
plot(t,yRK2,':')


%% RK2 stability considerations, FDE analysis
adt=linspace(0.01,3,20);
ladt=numel(adt);
G=zeros(ladt,1);
for igain=1:ladt
    G(igain)=(1-adt(igain)+1/2*adt(igain).^2);
end %for
figure;
plot(adt,G,'o')
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');


%% RK4 example; comparison against first and second order methods
yRK4=zeros(1,N);
yRK4(1)=y0;
for n=2:N
    dy1=dt*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt*fRK(t(n-1)+dt,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for
figure(1);
plot(t,yRK4,'^-')
legend('exact','fwd Eul','bwd Eul.','RK2','RK4')


%% RK2 and systems of equations, oscillating charge example
q=-1.6e-19;
m=9.1e-31;
B=10000e-9;
omega=q*B/m;
tmin=0;
tmax=2*2*pi/abs(omega);
t=linspace(tmin,tmax,50);
dt=t(2)-t(1);
lt=numel(t);

vx=zeros(1,lt);
vy=zeros(1,lt);
vy(1)=1e3;
vx(1)=1e3;
for n=2:lt
    %step x and y components together, this is the half update
    vxhalf=vx(n-1)+dt/2*(omega*vy(n-1));
    vyhalf=vy(n-1)-dt/2*(omega*vx(n-1));
    
    %now the full update
    vx(n)=vx(n-1)+dt*(omega*vyhalf);
    vy(n)=vy(n-1)-dt*(omega*vxhalf);    
end %for
figure;
ax=plotyy(t,vx,t,vy);
set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');


%integrate velocity to get position as a fn. of time
x=cumtrapz(t,vx);    %Matlab built-in for accumulating an integral value
y=cumtrapz(t,vy);
vz=1e3;
z=vz*t;

%comet plot demo
figure;
comet3(x,y,z)
set(gca,'FontSize',20);
xlabel('x');
ylabel('y');
zlabel('z');


%% Handling multiple time scales (ODE stiffness, book example from Gear's paper)
% dy/dx=-alpha*(y=F(t))+F'(t)
% y(t)=(y0-F(0))exp(-alpha*t)+F(t)
y0=1;
alpha=1000;
figure;

tsmin=0;
tsmax=4;
dts=0.002;
ts=tsmin:dts:tsmax;
ybar=(y0-2)*exp(-alpha*ts)+ts+2;
Ns=numel(ts);

plot(ts,ybar);
hold on;

yfwds=zeros(1,Ns);
yfwds(1)=y0;
ybwds=zeros(1,Ns);
ybwds(1)=y0;
for n=2:Ns
    yfwds(n)=yfwds(n-1)+dts*(-1000*(yfwds(n-1)-ts(n-1)-2)+1);
    ybwds(n)=(ybwds(n-1)+1000*ts(n-1)*dts+2001*dts)/(1+1000*dts);    
end %for
plot(ts,yfwds);
plot(ts,ybwds);
axis([0 4 1 6]);
set(gca,'FontSize',20);
legend('theory','fwd','bwd')


