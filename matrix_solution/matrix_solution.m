% A simple script to solve a system of equations, while checking various
% aspects of the matrix representation of that system


%% Illustrate the number of operations needed to implement Cramer's rule

n=1:10;
nops=(n-1).*factorial(n+1)+n;      %see book for justification
loglog(n,nops,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue');
xlabel('size of system (# of unknowns)');
ylabel('number of multiply/divides needed to solve');
title('Computational Cost of Cramer''s Rule')


%% Define a problem

% Example solved by hand in class
%   x1 +   x2 + 3*x3  = 2
% 5*x1 + 3*x2 +   x3 = 3
% 2*x1 + 3*x2 +   x3 = -1
A=[1, 1, 3; ...
   5, 3, 1; ...
   2, 3, 1];
b=[2; 3; -1];
x=A\b;

disp('x = ');
disp(x);


%% Illustrate vanilla forward elimination

n=6;              %system size
Bref=randn(n,n+1);    %augmented matrix containing RHS of system of equations, hopefully not singular since using randn...
%Bref=cat(2,A,b);
%n=size(A,1);
bref=Bref(:,n+1);   %RHS

%note that the elimination procedure modifies the matrix B
B=Bref;          %make a copy of the original since it will be change by this code
for ir1=2:n                                           %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:n                                     %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=B(ir2,ir1-1);                                    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        B(ir2,:)=B(ir2,:)-fact/B(ir1-1,ir1-1).*B(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements)
    end %for
end %for

disp('elim(B) = ');
disp(B);


%% Illustrate back substitution on B

n=size(B,1);                   %number of unknowns in the system
xsoln=zeros(n,1);              %space in which to store our solution vector
xsoln(n)=B(n,n+1)/B(n,n);      %finalized solution for last variable

%note that B is assumed to be upper triangular at this point
for ir1=n-1:-1:1              %iterate backwards from last equation using its value to solve others
    for ir2=ir1:-1:1          %this must also iterate backwards
        xsoln(ir2)=B(ir2,n+1);    %assume we're only dealing with a single right-hand side here.
        fact=B(ir2,ir2);      %diagonal element to be divided through doing subs for the ir2 row
        for ic=ir2+1:n
            xsoln(ir2)=xsoln(ir2)-B(ir2,ic)*xsoln(ic);
        end %for
        xsoln(ir2)=xsoln(ir2)/fact;     %divide once at the end to reduce number of ops
    end %for
end %for

disp('Elimination/back sub solution:  ');
disp(xsoln);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(Bref(1:n,1:n)\bref);


%% Check the Gaussian elimination function
[Bmod,ord]=elim(Bref(1:n,1:n),bref);

disp('Elimination with scaled pivoting:  ')
disp(Bmod(ord,:));

xgauss=backsub(Bmod(ord,:));

disp('Back sub solution using Gaussian elimination result:  ')
disp(xgauss);




