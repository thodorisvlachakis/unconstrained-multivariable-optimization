% This function is the implementation of Newton's Method while
% step g is defined by solving an inner optimization problem in every iteration which is to minimize the
% function f( x(k)+ g(k).*d(1), y(k) +g(k).*d(2) ) (we use functions with two
% variables, so I will imply hte method for two variables) by g(k) and this g(k) is the step we look for. We need, as inputs, 
% the function f we want to minimize, the starting point x0 ,the constant e
% which determines the termination of the algorithm.
% The function returns the final estimation for the minimum of f and the point of minimization. Also, the
% function returns a vector F in which we store the estimation for the
% minimum of f in every iteration and a vector X in which we store the
% estimation for the point of minimization in every iteration.

% !!!! Newton's Method uses the hessian of function f to compute the vector d in every
% iteration and it is needed the hessian of function f calculated at the
% searcing point be positive definite. This is not ensured and that's why Newton's Method may not reach out the desirable solution.  

function [MinimumValue, MinimizationPoint, F, X, Iterations] = NewtonMethodDescentStepByInnerOptimization(f,x0,e)

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

if(k>100)
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

eigenvalues=eig(D);
for i=1:length(eigenvalues)
    if(eigenvalues(i)<=0)
        if(k==1)
            message=sprintf('During %dst iteration the hessian of function f(x,y) calculated at the searching point is not positive definite. \n', k);
            fprintf(message);
            disp("Newton's Method violates the property of the iterative descent.");
            
            X=[x;y];
            F(k)=f(x(k),y(k));
            MinimumValue=vpa( f(x(k),y(k)) );
            MinimizationPoint=[x(k),y(k)];
            Iterations=k;

            return
        elseif (k==2)
            message=sprintf('During %dnd iteration the hessian of function f(x,y) calculated at the searching point is not positive definite. \n', k);
            fprintf(message);
            disp("Newton's Method violates the property of the iterative descent.");
            
            X=[x;y];
            F(k)=f(x(k),y(k));
            MinimumValue=vpa( f(x(k),y(k)) );
            MinimizationPoint=[x(k),y(k)];
            Iterations=k;

            return
        else
            message=sprintf('During %dth iteration the hessian of function f(x,y) calculated at the searching point is not positive definite. \n', k);
            fprintf(message);
            disp("Newton's Method violates the property of the iterative descent.");

            X=[x;y];
            F(k)=f(x(k),y(k));
            MinimumValue=vpa( f(x(k),y(k)) );
            MinimizationPoint=[x(k),y(k)];
            Iterations=k;

            return
        end
    end



d=-D*Gradientf(x(k),y(k));

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