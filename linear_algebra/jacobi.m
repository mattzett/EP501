function x=jacobi(x0,A,b,verbose)

%% Check the inputs
narginchk(3,4);
if nargin<4
    verbose=false;
end %if


%% Setup iterations
maxit=100;    %max number of iterations
conv=1e-6;    %relative change in residual for convergence
n=size(A,1);  %system size
residual=10*ones(n,1);
difftot=1e3+conv;   %max sure we enter iterations
x=x0;


%% Perform iterations
it=1;
while(difftot>conv && it<=maxit)
    difftotprev=difftot;
    resprev=residual;
    xprev=x;
    for i=1:n
        residual(i)=b(i);
        for j=1:n
            residual(i)=residual(i)-A(i,j)*xprev(j);
        end %for
        x(i)=xprev(i)+residual(i)/A(i,i);
    end %for
    difftot=sum(abs(residual-resprev));
    
    if (verbose)
        fprintf('x= ');
        for i=1:n
            fprintf('%f   ',x(i));
        end %for
        fprintf('\n');
        fprintf('it=%d; difftot = %e\n',it,difftot);
    end %if
    
    if (difftot>difftotprev)
        error('Solution appears to be diverging, check diagonal dominance...')
    end %if
    it=it+1;
end %while

if (it-1==maxit)
    warning('Solution may not have converged fully...')
end %if

end %function
