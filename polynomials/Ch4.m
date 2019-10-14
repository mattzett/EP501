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


%% Another take using Matlab built-in (direct fit)
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


%% Linear least squares example (line)
a=2;
b=3;
xdata=linspace(-5,5,25);
ydata=a+b*xdata+2*randn(size(xdata));
figure;
plot(xdata,ydata,'o','MarkerSize',20);

