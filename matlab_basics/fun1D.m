%% Objective function for 1D examples
function fval=fun1D(x,A,x0,sigmax)

  fval=A.*exp(-(x-x0).^2/2/sigmax);

end %function