% A script to demonstrate solutions to nonlinear equations
%
% requires:  objfun?.m


%% Common setup for closed domain problems
maxit=100;       %maximum number of iterations
f=@objfun1;      %set the function for which we are finding roots, change to illustrate different problems
minx=0;          %interval over which we are finding root (closed domain problems)
maxx=2*pi;
tol=1e-6;        %how close to zero we need to get to cease iterations
x=linspace(minx,maxx,24);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;


%% Plot function for which roots are to be found
figure(1);
plot(x,ygrid,'-');
xlabel('x');
ylabel('f(x)');
title('Objective function evaluated on grid')


%% Illustrate bracketing to find a single root on a given interval
it=1;
converged=false;    %make sure we enter iterations
a0=pi/4;
b0=pi;               %isolate one root for now
a=a0;
b=b0;
while (~converged && it<=maxit)
    it=it+1;
    
    c=(a+b)/2;
    aprev=a;           %save these for plotting
    bprev=b;
    if(f(a)*f(c)<0)    %c will be our new endpoint
        b=c;
        left=true;
    else               %a will be our new staring point
        a=c;
        left=false;
    end %if
    
    xnew=c;
    fval=f(xnew);
    converged=abs(fval)<tol;
        
    if (verbose)       %plot progress of the algorithm toward convergence
        figure(2);
        clf;
        hold on;
        plot(x,ygrid);
        ax=axis;
        axis([a0,b0,ax(3:4)]);
        if (left)
            plot(aprev,0,'k^','MarkerSize',10,'LineWidth',2);
            plot(bprev,0,'r^','MarkerSize',10,'LineWidth',2);
        else
            plot(aprev,0,'r^','MarkerSize',10,'LineWidth',2);
            plot(bprev,0,'k^','MarkerSize',10,'LineWidth',2);
        end %if
        plot(c,0,'ko','MarkerSize',10,'LineWidth',2);
        hold off;
        xlabel('x');
        ylabel('f(x)');
        title(sprintf('x = %f',xnew))
        pause;
    end %if
end %while
if (it==maxit)
    warning('Max number of iterations used...')
end %if
disp('Root value through bracketing:  ');
disp(xnew);
disp('Number of iterations required to reach tolerance:  ');
disp(it-1);


