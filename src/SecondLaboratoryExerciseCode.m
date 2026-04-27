% Second Laboratory Exercise

addpath(genpath(pwd));

% 1. I am going to plot the function f(x,y) I want to minimize. I will plot
% it for x∈ [-4, 4] and y∈ [-4, 4]. As we can see there is no reason to plot
% it for greater values of x and y. I tried it for x∈ [-400, 400] and y∈ [-400, 400]. We just decrease 
% the algorithm's performance. 

% Definition of the function f(x,y)
syms x y;
f(x,y)=(x.^3).* exp(-x.^2 -y.^4);

%fxy=@(x,y) (x.^3).* exp(-x.^2 -y.^4);

figure(1);
clf

fsurf(f,[-4 4 -4 4]);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Function f(x,y)');

% Before we start, we compute the exact minimum value of the function f by
% finding its critical points and compute the values at them.

fprintf('\n');

disp("Solving the differential equation ∇(f(x,y))=0 , the critical points of the function f are computed.");
Df(x,y)=gradient(f);
[X1, Y1]=solve( Df(x,y)==0, [x y]);

fprintf('The critical points are: ');
for i=1:length(X1)
    if(i>1)
        fprintf(', ');
    end

    fprintf('(%.10f , %.10f) ', X1(i), Y1(i));
end
fprintf('.\n');

disp("So these are the possible minimization points of the function f(x,y).");
fprintf('\n');
disp("Comparing the values of the function f at them, not only the minimization point is found, but also the minimum value of f(x,y).");

ExactPointOfMinimization=[X1(1) Y1(1)];
ExactMinimumValue= f( ExactPointOfMinimization(1), ExactPointOfMinimization(2) );

for i=1:length(X1)
    if( f( X1(i), Y1(i) )< ExactMinimumValue )
        ExactPointOfMinimization=[X1(i) Y1(i)];
        ExactMinimumValue= f( ExactPointOfMinimization(1), ExactPointOfMinimization(2) );

    end

end

fprintf('The exact minimization point of the function f(x,y) is: (%.10f , %.10f)\n', ExactPointOfMinimization(1), ExactPointOfMinimization(2) );
fprintf('The exact minimum value of the function f(x,y) is: %.10f\n',ExactMinimumValue);
fprintf('\n');
fprintf('\n');
% !! In the methods below the termination condition is |∇(f(x,y))|<e where
% e is small enough to satisfy the required accuracy on determination of minimization point
% and minimum value of the function f.
% So, firstly we will define that e small enough is order to take a good
% estimation while executing the algorithms below.
e=0.001;

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% 2. Steepest Descent Method
% We will execute the "Steepest Descent Method" algorithm starting from
% three different points (i.e. three difderent cases). In each case, the
% descent step g(k) will be chosen with three different ways, so there will be three different nested cases in each case.
fprintf('!!!! Steepest Descent Method !!!!\n');

% Firstly, we define the three starting points.
x01=[0 0];
x02=[-1 -1];
x03=[1 1];

% i) x01=[0 0]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Steepest Descent Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodWithConstantDescentStep (f,x01,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(2);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Steepest Descent Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodDescentStepByInnerOptimization (f,x01,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(3);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Steepest Descent Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodArmijoStepSizeRule (f,x01,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(4);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% ii) x02=[-1 -1]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Steepest Descent Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodWithConstantDescentStep (f,x02,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(5);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Steepest Descent Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodDescentStepByInnerOptimization (f,x02,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(6);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Steepest Descent Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodArmijoStepSizeRule (f,x02,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(7);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% iii) x03=[1 1]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Steepest Descent Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodWithConstantDescentStep (f,x03,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(8);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Steepest Descent Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodDescentStepByInnerOptimization (f,x03,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(9);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Steepest Descent Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Steepest Descent Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= SteepestDescentMethodArmijoStepSizeRule (f,x03,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(10);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% 3. Newton's Method
% We will execute the "Newton's Method" algorithm starting from
% three different points (i.e. three difderent cases). In each case, the
% descent step g(k) will be chosen with three different ways, so there will be three different nested cases in each case.
disp("!!!! Newton's Method !!!!");

% Firstly, we define the three starting points.
x01=[0 0];
x02=[-1 -1];
x03=[1 1];

% i) x01=[0 0]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Newton Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= NewtonMethodWithConstantDescentStep (f,x01,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(11);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Newton Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= NewtonMethodDescentStepByInnerOptimization (f,x01,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(12);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Newton Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= NewtonMethodArmijoStepSizeRule (f,x01,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(13);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% ii) x02=[-1 -1]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Newton Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, ~, ~, Iterations]= NewtonMethodWithConstantDescentStep (f,x02,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.
% Because of non-convergence of the function f we use the amendment in
% order to prove it.
[~, ~, F, ~, Iterations]= NewtonMethodWithConstantDescentStepAmendment (f,x02,e,g0);

figure(14);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Newton Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, ~, ~, Iterations]= NewtonMethodDescentStepByInnerOptimization (f,x02,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.
% Because of non-convergence of the function f we use the amendment in
% order to prove it.
[~, ~, F, ~, Iterations]= NewtonMethodWithConstantDescentStepAmendment (f,x02,e,g0);

figure(15);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Newton Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, ~, ~, Iterations]= NewtonMethodArmijoStepSizeRule (f,x02,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.
% Because of non-convergence of the function f we use the amendment in
% order to prove it.
[~, ~, F, ~, Iterations]= NewtonMethodArmijoStepSizeRuleAmendment (f,x02,e,s,a,b);

figure(16);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% iii) x03=[1 1]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Newton Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, ~, ~, Iterations]= NewtonMethodWithConstantDescentStep (f,x03,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.
% Because of non-convergence of the function f we use the amendment in
% order to prove it.
[~, ~, F, ~, Iterations]= NewtonMethodWithConstantDescentStepAmendment (f,x03,e,g0);

figure(17);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Newton Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, ~, ~, Iterations]= NewtonMethodDescentStepByInnerOptimization (f,x03,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.
% Because of non-convergence of the function f we use the amendment in
% order to prove it.
[~, ~, F, ~, Iterations]= NewtonMethodDescentStepByInnerOptimizationAmendment (f,x03,e);

figure(18);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Newton Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Newton Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, ~, ~, Iterations]= NewtonMethodArmijoStepSizeRule (f,x03,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.
% Because of non-convergence of the function f we use the amendment in
% order to prove it.
[~, ~, F, ~, Iterations]= NewtonMethodArmijoStepSizeRuleAmendment (f,x03,e,s,a,b);

figure(19);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% 4. Levenberg-Marquadt Method
% We will execute the "Levenberg-Marquadt Method" algorithm starting from
% three different points (i.e. three difderent cases). In each case, the
% descent step g(k) will be chosen with three different ways, so there will be three different nested cases in each case.
disp("!!!! Levenberg-Marquadt Method !!!!");

% Firstly, we define the three starting points.
x01=[0 0];
x02=[-1 -1];
x03=[1 1];

% i) x01=[0 0]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Levenberg-Marquadt Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, X, Iterations]= LevenbergMarquadtMethodWithConstantDescentStep (f,x01,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(20);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Levenberg-Marquadt Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodDescentStepByInnerOptimization (f,x01,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(21);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x01(1), x01(2));
fprintf('Execution of Levenberg-Marquadt Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodArmijoStepSizeRule (f,x01,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(22);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% ii) x02=[-1 -1]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Levenberg-Marquadt Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodWithConstantDescentStep (f,x02,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(23);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Levenberg-Marquadt Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodDescentStepByInnerOptimization (f,x02,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(24);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x02(1), x02(2));
fprintf('Execution of Levenberg-Marquadt Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodArmijoStepSizeRule (f,x02,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(25);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% iii) x03=[1 1]. We will execute the algorithm and we will take the results
% and the graph of convergence of the function f as a function of iterations.

% a) In the first subcase the descent step g(k) will be constant for every
% k (iteration). We choose g(k)=g0=0.1 which is a value for g(k) that
% ensures the convergence of the method, according to Theorem 5.2.6 in
% book.

g0=0.1;

fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Levenberg-Marquadt Method with constant descent step g(k)=%f for every k.\n',g0);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodWithConstantDescentStep (f,x03,e,g0);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(26);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');

% b) In the second subcase the descent step g(k) will be determined by solving an inner optimization problem
% in every iteration which is to minimize the function f( x(k)+ g(k).*d(1),
% y(k) +g(k).*d(2) ) by g(k) and this g(k) is the step we look for.
% (d=∇f(x(k),y(k)) for the kth iteration of Steepest Descent Method). 
fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Levenberg-Marquadt Method while descent step g(k) will be determined by solving an inner optimization problem in every iteration which is to minimize the function f(x(k+1)) by g(k).\n');
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodDescentStepByInnerOptimization (f,x03,e);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(27);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');


% c) In the third subcase the descent step g(k) will be determined by Armijo Step Size Rule in every iteration as we have implied it.
% Armijo Step Size Rule needs parameters a,b,g0 to be defined. We will
% choose a=0.01 , b=0.3, g0=0.25 which are well defined in order to run the
% algorithm while step g(k) is defined by Armijo Step Size Rule in every
% iteration.

s=0.25;
% Use s instead of g0.
a=0.01;
b=0.3;

fprintf('\n');
fprintf('Levenberg-Marquadt Method: The starting point is the point (%f , %f) .\n',x03(1), x03(2));
fprintf('Execution of Levenberg-Marquadt Method while descent step g(k) will be determined by Armijo Step Size Rule in every iteration.\n');
fprintf('Armijo Step Size Rule needs parameters a,b,s to be defined. So a well defined choice is a=%f , b=%f , s=%f .\n',a,b,s);
fprintf('\n');

[MinimumValue, MinimizationPoint, F, ~, Iterations]= LevenbergMarquadtMethodArmijoStepSizeRule (f,x03,e,s,a,b);

fprintf('The algorithm is executed for %d Iteration(s) ending up at the point (%.10f , %.10f) .\n',Iterations, MinimizationPoint(1), MinimizationPoint(2));
fprintf('So, the estimation for the Minimization Point is: (x*,y*)=(%.10f , %.10f) and for the Minimum Value of the function f is: f(x*,y*)=%.10f .\n',MinimizationPoint(1), MinimizationPoint(2), MinimumValue);

% Graph of convergence of the function f as a function of iterations.

figure(28);
clf

k=1:Iterations;
scatter(k, F,30, "blue","filled");
xlabel('Iterations');
ylabel('f(x(k),y(k))');
title('Graph of convergence of function f(x,y) as a function of Iterations');
