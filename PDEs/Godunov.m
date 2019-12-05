function fnew=Godunov(dt,dx,v,f)

% Implements directional-based upwinding (Godunov's method) for solving
% hyperbolic PDES.

lx=numel(f);
fleft=f(lx-1:lx);    %ghost cells here implement periodic boundary conditions
fright=f(1:2);

fnew=zeros(size(f));
if (v<0)   %upwind direction is i+1, since travelling backwards
    fnew(1)=f(1)-v*dt/dx*(f(2)-f(1));
else
    fnew(1)=f(1)-v*dt/dx*(f(1)-fleft(2));
end %if
for ix=2:lx-1
    if (v<0)
        fnew(ix)=f(ix)-v*dt/dx*(f(ix+1)-f(ix));
    else
        fnew(ix)=f(ix)-v*dt/dx*(f(ix)-f(ix-1));
    end %if
end %for
if (v<0)   %upwind direction is i+1, since travelling backwards
    fnew(lx)=f(lx)-v*dt/dx*(fright(1)-f(lx));
else
    fnew(lx)=f(lx)-v*dt/dx*(f(lx)-f(lx-1));
end %if

end %function