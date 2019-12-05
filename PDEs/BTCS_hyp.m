function fnew=BTCS_hyp(dt,dx,v,f)

% Implements the BTCS method for solving hyperbolic equations.
% Performs a single time update for time step dt.  By default this code
% will assume periodic boundary conditions.  

lx=numel(f);
fleft=f(lx-1:lx);    %ghost cells here implement periodic boundary conditions
fright=f(1:2);

fnew=zeros(size(f));

%Define our system of equations
A=sparse(lx,lx);
b=zeros(lx,1);

%leftmost grid point
A(1,lx)=-v/2/dx;
A(1,1)=1/dt;
A(1,2)=v/2/dx;
b(1)=f(1)/dt;

%interior grid points
for ix=2:lx-1
    %i-1 term in i equation
    A(ix,ix-1)=-v/2/dx;
    
    %i
    A(ix,ix)=1/dt;
    
    %i+1
    A(ix,ix+1)=v/2/dx;
    
    %RHS
    b(ix)=f(ix)/dt;
end %for

%rightmost grid point, this creates a circulant matrix, which matlab seems
%perfectly happy to deal with...
A(lx,lx-1)=-v/2/dx;
A(lx,lx)=1/dt;
A(lx,1)=v/2/dx;
b(lx)=f(lx)/dt;

%Matlab sparse matrix solution
fnew=A\b;

end %function