% A simple script to solve a system of equations, while checking various
% aspects of the matrix representation of that system
%
%This script requires files containing the functions Gauss_elim, backsub, and Jacobi.  


%% Illustrate the number of operations needed to implement Cramer's rule
n=1:10;
nops=(n-1).*factorial(n+1)+n;      %see book for justification
figure(1);
loglog(n,nops,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue');
xlabel('size of system (# of unknowns)');
ylabel('number of multiply/divides needed to solve');
title('Theoretical Computational Cost of Cramer''s Rule')


%% Define the example problem discussed in class
% Example solved by hand in class
%   x1 + 4*x2 + 2*x3  = 15
% 3*x1 + 2*x2 +   x3 = 10
% 2*x1 +   x2 + 3*x3 = 13
A=[1, 4, 2; ...
   3, 2, 1; ...
   2, 1, 3];
b=[15;10;13];
x=A\b;
disp('(class problem Matlab solution) x = ');
disp(x);


%% Illustrate vanilla forward elimination
nref=6;                %system size for larger reference problem
Aref=randn(nref,nref);    %augmented matrix containing RHS of system of equations, in practice you'd want to check conditioning...
bref=randn(nref,1);    %RHS

%note that the elimination procedure coded below modifies the matrix B
Awork=cat(2,Aref,bref);          %This is our working version of the matrix used to perform elimination (i.e. it will be modified)
for ir1=2:nref                                           %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:nref                                     %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=Awork(ir2,ir1-1);                                    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        Awork(ir2,:)=Awork(ir2,:)-fact/Awork(ir1-1,ir1-1).*Awork(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements), this is a little bit wasteful as it uses entire row...
    end %for
end %for

disp('elim([Aref,bref]) = ');
disp(Awork);


%% Illustrate back substitution on B using provided Matlab function
xsoln=backsub(Awork);
disp('Elimination/back sub solution:  ');
disp(xsoln);
disp('Matlab,GNU/Octave built-in solution:  ');
disp(Aref\bref);


%% Use the Gaussian elimination function to solve the same system (include scaled pivoting)
[Amod,ord]=Gauss_elim(Aref,bref);

disp('Elimination with scaled pivoting on matrix:  ');
disp(Amod(ord,:));
xgauss=backsub(Amod(ord,:));
disp('Back substitution solution using Gaussian elimination result:  ');
disp(xgauss);


%% Print step by step solution (Gauss elimination) for a simple system to illustrate
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
[Amodsmall,ord]=Gauss_elim(A,b,true);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


%% Solve a diagonally dominant system using Jacobi iteration
nit=10;
Ait=diag(-1*ones(nit-1,1),-1)+diag(-1*ones(nit-1,1),1)+diag(4*ones(nit,1),0);    %this must be diagonally dominant or else the method won't converge
%Ait=randn(nit,nit);    %see if code can detect non-diagonal dominance and exit gracefully...
x0=randn(nit,1);
bit=ones(nit,1);
tol=1e-9;
disp('Verbose Jacobi iterations:  ')
[xit,iterations]=Jacobi(x0,Ait,bit,tol,true);

disp('Solution with Jacobi iterations:  ')
disp(xit);
disp('Number of iterations required and tolerance:  ')
disp(iterations);
disp(tol);
disp('Matlab built-in solution:  ')
disp(Ait\bit);


%% Evaluate performance and scaling of Gaussian elimination by solving systems of different size
nvals=50:50:500;
testtimes=zeros(size(nvals));
lrep=1;     %how many times to repeat each test

disp('Start of tests of Gaussian-elimination scaling');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=randn(nlarge,nlarge);
    blarge=randn(nlarge,1);
    
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        [Blargemod,ordlarge]=Gauss_elim(Blarge,blarge);
        xlarge=backsub(Blargemod(ordlarge,:));
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(2);
plot(nvals,testtimes,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
title('Empirical Performance of Gaussian Elimination');
