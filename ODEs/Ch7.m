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
figure;
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);


%% Forward Euler solution
yfwd=zeros(1,N);
yfwd(1)=y0;
for it=2:N
    yfwd(it)=yfwd(it-1)*(1-alpha*dt);
end %for
hold on;
plot(t,yfwd,'--');


%% Backward Euler solution
ybwd=zeros(1,N);
ybwd(1)=y0;
for it=2:N
    ybwd(it)=ybwd(it-1)/(1+alpha*dt);
end %for
plot(t,ybwd,'-.');
legend('exact','fwd Eul','bwd Eul.')


%% Illustration of stability
alpha2=10;
dtmargin=2/alpha2;
dtstable=1/10*dtmargin;
dtunstable=1.01*dtmargin;

dt=dtunstable;
t2=0:dt:5;
N=numel(t);
y2bar=y0*exp(-alpha2*t2);

y2fwd=zeros(1,N);
y2fwd(1)=y0;
for it=2:N
    y2fwd(it)=y2fwd(it-1)*(1-alpha2*dt);
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


%% RK2 stability considerations, FDE analysis


%% RK4 example; comparison against first and second order methods


%% RK2 and systems of equations, oscillating charge example


%% Handling multiple time scales (ODE stiffness, book examples)



