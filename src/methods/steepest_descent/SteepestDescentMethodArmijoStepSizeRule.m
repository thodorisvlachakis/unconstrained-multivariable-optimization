% This function is the implementation of Steepest Descent Method while
% step g is defined by Armijo Step Size Rule in every iteration. The step g0 is the initial
% step we use in Armijo function and this has to be given as an input of the function. We also need, as inputs, 
% the function f we want to minimize, the starting point x0 ,the constant e
% which determines the termination of the algorithm and the parametres a and b which are used in Armijo function.
% The function returns the final estimation for the minimum of f and the point of minimization. Also, the
% function returns a vector F in which we store the estimation for the
% minimum of f in every iteration and a vector X in which we store the
% estimation for the point of minimization in every iteration.


function [MinimumValue, MinimizationPoint, F, X, Iterations] = SteepestDescentMethodArmijoStepSizeRule(f,x0,e,g0,a,b)

k=1;

% I will define two vectors x,y in which i will store the coordinats x and
% y,respectively, in every iteration. 
x(k)=x0(1);
y(k)=x0(2);

% I need an array F in which I will store the values of function f in every
% iteration.
% Also, I create two vectors m and g where the elements m(k) and g(k) ,respectively, are the outputs of the
% AmrijoStepSizeRule function which is used into while loop, for the kth
% iteration.
F=zeros;
m=zeros;
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

[DescentStep, IntegerM]=ArmijoStepSizeRule(f,g0,a,b,d,x,y,k);
g(k)=DescentStep(1);
m(k)=IntegerM(1);

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