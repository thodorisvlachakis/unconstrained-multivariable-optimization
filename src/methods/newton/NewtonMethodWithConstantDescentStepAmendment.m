% Amendment of Newton's Algorithm in order to prove that this method is not
% converging at the minimization point of function f we want to minimize.
% This function refers to constant descent step g(k)=g0.

function [MinimumValue, MinimizationPoint, F, X, Iterations]= NewtonMethodWithConstantDescentStepAmendment (f,x0,e,g)

k=1;

% I will define two vectors x,y in which i will store the coordinats x and
% y,respectively, in every iteration. 
x(k)=x0(1);
y(k)=x0(2);

% I need an array F in which I will store the values of function f in every
% iteration.
F=zeros;

Gradientf=gradient(f);
NormGradientf=norm(Gradientf);

if(NormGradientf(x(k),y(k))==0)
    % This means that the starting point is critical point of the function
    X=[x;y];
    F(k)=f(x(k),y(k));
    MinimumValue=vpa( f(x(k),y(k)) );
    MinimizationPoint=[x(k),y(k)];
    Iterations=k;

    return;
end

while( NormGradientf(x(k),y(k)) >=e )

if(k>40)
    X=[x;y];
    F(k)=f(x(k),y(k));
    MinimumValue=vpa( f(x(k),y(k)) );
    MinimizationPoint=[x(k),y(k)];
    Iterations=k;

    disp("The method is not efficient because of too many iterations.");
    return

end

point=[x(k) y(k)];    
D=inv( HessianOfFunctionAtSpecificPoint(f,point) );
% it is problable the matrix D not be positive definite so Newton's Method violates
% the property of the iterative descent.


d=-D*Gradientf(x(k),y(k));

x(k+1)=x(k)+(g(1).*d(1));
y(k+1)=y(k)+(g(1).*d(2));

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