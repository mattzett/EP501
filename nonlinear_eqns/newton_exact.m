function [root,it]=newton_exact(f,fprime,x0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f (input as a handle)
% given a function which computes the derivative

%% Error checking of input
narginchk(3,6);   %check for correct number of inputs to function
if (nargin<4)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<5)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<6)
    verbose=false;
end %if


%% Newton iterations
it=1;
root=x0;
fval=100*tol;
while(abs(fval)>tol && it<=maxit)
    root=root-f(root)./fprime(root);    % update root estimate
    fval=f(root);                          % see how far off we are from zero...
    if (verbose)
        fprintf('iteration: %d; root:  %f; function value: %f \n',it,root,fval);
    end %if
    it=it+1;
end %while
it=it-1;

if (it==maxit)
    warning('Used max number of iterations, results are probably not good...')
end %if

end %function