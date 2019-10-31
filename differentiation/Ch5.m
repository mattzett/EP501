%% Some examples of numerical differentiation
lx=25;
x=linspace(-10,10,lx)';
y=sin(0.5*x);
yprime=0.5*cos(0.5*x);

figure;
plot(x,y);
hold on;
plot(x,yprime);


%% Comparison of basic numerical derivative
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


%first order derivative approximation (backward)
%interior
dy_dxbwd=zeros(lx,1);
for ix=2:lx
    dy_dxbwd(ix)=(y(ix)-y(ix-1))/dx;
end %for
dy_dxbwd(1)=dy_dxbwd(2);

plot(x,dy_dxbwd,'m--')
legend('original function','analytical','centered','backward')
xlabel('x');
ylabel('y(x) or y''(x)');
title('Comparison of finite difference derivative approximations');


%% Demonstrate effects of grid refinement (convergence)
lx2=256;
x2=linspace(-10,10,lx2)';
y2=sin(0.5*x2);
y2prime=0.5*cos(0.5*x2);


%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)

