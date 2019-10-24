%% Some examples of numerical differentiation
lx=25;
x=linspace(-10,10,lx)';
y=sin(0.5*x);
yprime=0.5*cos(0.5*x);

figure;
plot(x,y);
hold on;
plot(x,yprime);


%% Compute a numerical derivative
%second order, centered
dy_dx=zeros(lx,1);
dx=x(2)-x(1);

%forward difference at the beginning
dy_dx(1)=(y(2)-y(1))/dx;

%centered difference on the interior
for ix=2:lx-1
    dy_dx(ix)=(y(ix+1)-y(ix-1))/2/dx;
end %for

%backward difference at the end
dy_dx(lx)=(y(lx)-y(lx-1))/dx;

plot(x,dy_dx,'k--')


%first order (backward)
%backward difference 

%interior
for ix=2:lx
    dy_dx(ix)=(y(ix)-y(ix-1))/dx;
end %for
dy_dx(1)=dy_dx(2);

plot(x,dy_dx,'m--')


%% Some examples of numerical integration
