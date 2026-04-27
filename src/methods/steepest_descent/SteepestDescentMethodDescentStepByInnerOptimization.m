% This function is the implementation of Steepest Descent Method while
% step g is defined by solving an inner optimization problem in every iteration which is to minimize the
% function f( x(k)+ g(k).*d(1), y(k) +g(k).*d(2) ) (we use functions with two
% variables, so I will imply hte method for two variables) by g(k) and this g(k) is the step we look for. We need, as inputs, 
% the function f we want to minimize, the starting point x0 ,the constant e
% which determines the termination of the algorithm.
% The function returns the final estimation for the minimum of f and the point of minimization. Also, the
% function returns a vector F in which we store the estimation for the
% minimum of f in every iteration and a vector X in which we store the
% estimation for the point of minimization in every iteration.


function [MinimumValue, MinimizationPoint, F, X, Iterations] = SteepestDescentMethodDescentStepByInnerOptimization(f,x0,e)

k=1;

% I will define two vectors x,y in which i will store the coordinats x and
% y,respectively, in every iteration. 
x(k)=x0(1);
y(k)=x0(2);

% I need an array F in which I will store the values of function f in every
% iteration.

F=zeros;
g=zeros;

Gradientf=gradient(f);
NormGradientf=norm(Gradientf);

while( NormGradientf(x(k),y(k)) >=e )
   
if(k>100)
    X=[x;y];
    F(k)=f(x(k),y(k));
    MinimumValue=vpa( f(x(k),y(k)) );
    MinimizationPoint=[x(k),y(k)];
    Iterations=k;

    disp("The method is not efficient because of too many iterations.");
    return

end    
    
    
d=-Gradientf(x(k),y(k));

DescentStep = DescentStepByInnerOptimization(f,x,y,k,d);
g(k)=DescentStep(1);

x(k+1)=x(k)+(g(k).*d(1));
y(k+1)=y(k)+(g(k).*d(2));

F(k)=f(x(k),y(k));

k=k+1;

end


F(k)=f(x(k),y(k));

Iterations=k;


% I return the vectors x and y creating a 2×k matrix which has x as first
% row an y as second.
X=[x;y];

% I return the final estimations. 
MinimumValue=vpa( f(x(k),y(k)) );

MinimizationPoint=[x(k),y(k)];


end