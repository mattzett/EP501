%% Generate data from known polynomial and direct fit it for coeffs., verify
x=[1,2,3,4]';
y=2.*x.^3-3.*x.^2+4.*x+9;
figure;
plot(x,y,'*','MarkerSize',20);


%% Matlab built-in (direct fit)
n=3;
coeffs=polyfit(x,y,n);
xlarge=linspace(1,4,50);
ylarge=polyval(coeffs,xlarge);
hold on;
plot(xlarge,ylarge,'--');


%% Another take using Matlab built-in (a direct fit)
X=zeros(n+1,n+1);
for icol=n+1:-1:1
    newcol=x(:).^(icol-1);
    X(:,n+2-icol)=newcol;
end %for
coeffs2=X\y;

ylarge2=polyval(coeffs2,xlarge);
plot(xlarge,ylarge2,'.','MarkerSize',20);


%% Can use repo Gauss elimination to solve (direct fit)
addpath ../linear_algebra/;
[Xmod,ord]=Gauss_elim(X,y);
coeffs3=backsub(Xmod(ord,:));
rmpath ../linear_algebra/; 


%% Quick demo of Matlab random number generation features
figure;
uniform_noise=rand([10000,1]);
hist(uniform_noise,25);
xlabel('random number');
ylabel('# of occurences');

figure;
gaussian_noise=randn([10000,1]);
hist(gaussian_noise,25);
xlabel('random number');
ylabel('# of occurences');

%Create noise with a given standard deviation and mean
xmean=5;
xstdev=3;
gaussian_noise2=xstdev.*gaussian_noise+xmean;

figure;
hist(gaussian_noise2,25);
xlabel('random number');
ylabel('# of occurences');


%% Linear least squares example (fitting to a line - a first-order poly)
% y=a+b*x
n=1000;
a=2;
b=3;
xdata=linspace(-5,5,n);         %independent variable
ynoise=8*randn(size(xdata));     %noise with standard dev. of 2, zero mean.
ytrue=a+b*xdata;
ydata=ytrue+ynoise;

%Plot of noisy data
figure;
plot(xdata,ytrue,'k-');
xlabel('x');
ylabel('y');
title('Illustrating a linear fit')
hold on;
plot(xdata,ydata,'o','MarkerSize',20);

%Set up the Jacobian for an elimination fit to a line
addpath ../linear_algebra/;
 
J=cat(2,ones(n,1),xdata(:));
M=J'*J;
yprime=J'*ydata(:);
[Mmod,ord]=Gauss_elim(M,yprime);
avec=backsub(Mmod(ord,:));
yfit=avec(1)+avec(2)*xdata;

rmpath ../linear_algebra/;

plot(xdata,yfit,'--');
legend('true function','data','linear fit');
hold off;


%% Illustration of bilinear interpolation
