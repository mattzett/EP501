function x=backsub(A)

% This function performs back substitution on an upper traingular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be upper triangular at this point.


n=size(A,1);                   %number of unknowns in the system
x=zeros(n,1);              %space in which to store our solution vector
x(n)=A(n,n+1)/A(n,n);      %finalized solution for last variable, resulting for upper triangular conversion

for ir1=n-1:-1:1              %iterate backwards from last equation using its value to solve others
    for ir2=ir1:-1:1          %this must also iterate backwards
        x(ir2)=A(ir2,n+1);    %assume we're only dealing with a single right-hand side here.
        fact=A(ir2,ir2);      %diagonal element to be divided through doing subs for the ir2 row
        for ic=ir2+1:n
            x(ir2)=x(ir2)-A(ir2,ic)*x(ic);
        end %for
        x(ir2)=x(ir2)/fact;     %divide once at the end to reduce number of ops
    end %for
end %for

end %function