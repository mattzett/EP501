%% Objective function for 2D examples
function fval=fun2D(x,y,A,x0,y0,sigmax,sigmay)

  fval=A.*exp(-(x-x0).^2/2/sigmax^2).*exp(-(y-y0).^2/2/sigmay^2);

end %function