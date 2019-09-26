% A script to demonstrate solutions to nonlinear equations on closed
% domains
%
% requires:  objfun?.m (set function pointer f to desired function at beginning of program)


%% Common setup for closed domain problems
maxit=100;       %maximum number of iterations allowed
f=@objfun2;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun2_deriv;
minx=0;          %interval over which we are finding root (closed domain problems)
maxx=2*pi;
tol=1e-9;        %how close to zero we need to get to cease iterations
x=linspace(minx,maxx,24);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;


%% Plot function for which roots are to be found
figure(1);
plot(x,ygrid,'-');
xlabel('x');
ylabel('f(x)');
title('Objective function evaluated on grid')


%% Illustrate inverval halving to find a single root on a given interval
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
disp('Root value through inverval halving:  ');
disp(xnew);
disp('Number of iterations required to reach tolerance:  ');
disp(it-1);


%% Illustrate false position to find a single root on a given interval
it=1;
converged=false;    %make sure we enter iterations
a=a0;
b=b0;
while (~converged && it<=maxit)
    it=it+1;
    
    c=a-f(a)/((f(b)-f(a))/(b-a));
    
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
        y=(f(bprev)-f(aprev))/(bprev-aprev)*(x-aprev)+f(aprev);    %approximate line
        plot(x,y,'--');
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
disp('Root value through false position:  ');
disp(xnew);
disp('Number of iterations required to reach tolerance:  ');
disp(it-1);


%% Newton-Rhapson root-finding method
[xNewton,itNew]=newton_exact(f,fprime,0,100,tol,verbose);
disp('Root value through Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);

[xNewton,itNew]=newton_exact(f,fprime,4,100,tol,verbose);
disp('Root value through Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);


%% Newton-Rhapson for finding *all* roots



