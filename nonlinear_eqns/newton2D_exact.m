function [rootx,rooty,it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f,g (input as a handle)
% given a function which computes the derivative
% 
% requires:  Gauss_elim.m and backsub.m


%% Need access to linear algebra routines for solves at each iteration
addpath ../linear_algebra;


%% Error checking of input
narginchk(6,9);   %check for correct number of inputs to function
if (nargin<7)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<8)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<9)
    verbose=false;
end %if


%% Make sure we don't start at an inflection point with zero derivative
[gradfx,gradfy]=gradf(x0,y0);
[gradgx,gradgy]=gradg(x0,y0);
if (abs(min([gradfx,gradfy,gradgx,gradgy]))<tol)    %this needs to really check inflection vs. saddle point; will fix later
    warning(' Attempting to start Newton iterations near an inflection/saddle point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
    y0=y0+1;
end %if


%% Newton iterations
it=1;
rootx=x0;
rooty=y0;
fval=f(rootx,rooty);
gval=g(rootx,rooty);
converged=false;
while(~converged && it<=maxit)
    [gradfx,gradfy]=gradf(rootx,rooty);
    [gradgx,gradgy]=gradg(rootx,rooty);
    A=[gradfx,gradfy; ...
        gradgx,gradgy];
    
    fvec=[fval;gval];
    [Amod,ord]=Gauss_elim(A,-1*fvec);
    dxvec=backsub(Amod(ord,:));
    detA=prod(diag(Amod(ord,:)));
    if (abs(detA) < 1e-6)
        error(' Ended up at a point where Newton iteration is singular, try a different starting point')
    end %if
    
    rootx=rootx+dxvec(1);
    rooty=rooty+dxvec(2);
    fval=f(rootx,rooty);                  % see how far off we are from zero...
    gval=g(rootx,rooty);
    if (verbose)
        fprintf(' iteration: %d; x:  %f + %f i; y:  %f + %f i; f: %f, g:  %f\n',it,real(rootx),imag(rootx),real(rooty),imag(rooty),fval,gval);
    end %if
    it=it+1;
    converged=abs(fval)+abs(gval)<tol;
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function