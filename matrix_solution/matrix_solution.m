% A simple script to solve a system of equations, while checking various
% aspects of the matrix representation of that system


%% Illustrate the number of operations needed to implement Cramer's rule

n=[1,3,10,30,100,300,1000];
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
B=randn(n,n+1);    %augmented matrix containing RHS of system of equations, hopefully not singular since using randn...

for ir1=2:n                                           %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from for this particular column
    for ir2=ir1:n                                     %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=B(ir2,ir1-1);                                    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        B(ir2,:)=B(ir2,:)-fact/B(ir1-1,ir1-1).*B(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements)
    end %for
end %for

disp('elim(B) = ');
disp(B);

