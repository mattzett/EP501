%% Some examples of numerical differentiation
lx=20;
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
set(gca,'FontSize',24);

% %% Demonstrate effects of grid refinement (convergence)
% lx2=256;
% x2=linspace(-10,10,lx2)';
% y2=sin(0.5*x2);
% y2prime=0.5*cos(0.5*x2);


%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)
lx=20;
ly=20;
x=linspace(-5,5,lx);
y=linspace(-5,5,ly);
[X,Y]=meshgrid(x,y);
f=exp(-(X.^2)/2/2).*exp(-Y.^2/2/1);

figure;
contourf(x,y,f);

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);


%% Gradient
gradx=zeros(size(f));
grady=zeros(size(f));

%x component of gradient
for ix=2:lx-1
    gradx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;    %\partial/\partial x
end %for
gradx(:,1)=(f(:,2)-f(:,1))/dx;
gradx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;

%y component of gradient
for iy=2:ly-1
    grady(iy,:)=(f(iy+1,:)-f(iy-1,:))/2/dy;    %\partial/\partial y
end %for
grady(1,:)=(f(2,:)-f(1,:))/dy;
grady(ly,:)=(f(ly,:)-f(ly-1,:))/dy;

%add quiver on top of color plot
hold on;
quiver(X,Y,gradx,grady,'Color','white','LineWidth',2);
set(gca,'FontSize',24);


%% Take the Laplacian by taking divergence of the previously computed gradient
f=gradx;
g=grady;

%x derivative part of the divergence
divx=zeros(size(f));
for ix=2:lx-1
    divx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;
end %for
divx(:,1)=(f(:,2)-f(:,1))/dx;
divx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;

%y derivative part of the divergence
divy=zeros(size(y));
for iy=2:ly-1
    divy(iy,:)=(g(iy+1,:)-g(iy-1,:))/2/dy;
end %for
divy(1,:)=(g(2,:)-g(1,:))/dy;
divy(ly,:)=(g(ly,:)-g(ly-1,:))/dy;

div=divx+divy;    %this is really laplacian b/c input is gradient

figure;
surface(x,y,div);
set(gca,'FontSize',24);

