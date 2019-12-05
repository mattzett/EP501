function fnew=LaxFried(dt,dx,v,f)

% Implements the Lax-Friedrich's method for solving hyperbolic equations.
% Performs a single time update for time step dt.  By default this code
% will assume periodic boundary conditions.  

lx=numel(f);
fleft=f(lx-1:lx);    %ghost cells here implement periodic boundary conditions
fright=f(1:2);

fnew=zeros(size(f));
fnew(1)=1/2*(fleft(2)+f(2))-dt/2/dx*v*(f(2)-fleft(2));
for i=2:lx-1     %interior grid points
    fnew(i)=1/2*(f(i-1)+f(i+1))-dt/2/dx*v*(f(i+1)-f(i-1));
end %for
fnew(lx)=1/2*(f(lx-1)+fright(1))-dt/2/dx*v*(fright(1)-f(lx-1));

end %function