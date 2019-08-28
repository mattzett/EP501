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

