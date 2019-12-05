function fnew=LaxWen(dt,dx,v,f)

% Implements the Lax-Wendroff method for solving hyperbolic equations.
% Performs a single time update for time step dt.  By default this code
% will assume periodic boundary conditions. This particular code implements
% the two-step approach to LW.

lx=numel(f);
fleft=f(lx-1:lx);    %ghost cells here implement periodic boundary conditions
fright=f(1:2);

%half step lax-f update for cell edges. note indexing for fhalf,
%i->i-1/2, i+1->i+1/2
fhalf(1)=1/2*(fleft(2)+f(1))-dt/2/dx*v*(f(1)-fleft(2));
for i=2:lx     %interior grid points
    fhalf(i)=1/2*(f(i-1)+f(i))-dt/2/dx*v*(f(i)-f(i-1));
end %for
fhalf(lx+1)=1/2*(f(lx)+fright(1))-dt/2/dx*v*(fright(1)-f(lx));

%full time step LW update
fnew=zeros(size(f));
fnew(1)=f(1)-dt/dx*v*(fhalf(2)-fhalf(1));
for i=2:lx-1     %interior grid points
    fnew(i)=f(i)-dt/dx*v*(fhalf(i+1)-fhalf(i));
end %for
fnew(lx)=f(lx)-dt/dx*v*(fhalf(lx+1)-fhalf(lx));

end %function