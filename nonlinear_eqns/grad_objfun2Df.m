function [fx,fy]=grad_objfun2Df(x,y)

fx=3.*x.^2-3.*y;
fy=3.*y.^2-3.*x;

end %function