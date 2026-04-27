% This function is an implementation of inner optimization which can be used on each
% iteration of an algorithm in order to find the descent step g(k) (which regards to kth iteration of main algorithm)
% we need to use, so that we determine the next searching point x(k+1) and
% y(k+1).
% This step g(k) is computed by solving an inner optimization problem: The
% step g(k) is the one that minimize the function f( x(k)+ g(k).*d(1), y(k) +g(k).*d(2) ) (we use functions with two
% variables, so I will imply the method for two variables) which is a
% function of g(k). So, on kth iteration of main algorithm, we have to
% solve an inner optimization problem which is a problem of minimization of
% a function with only one variable. This problem can be solved by methods for minimization of such a function
% (e.g Bisection Method Using Derivatives or Golden Section Method). I am
% going to use Golden Section Method as the method of inner
% optimization. 
% There is a theorem which can be used to prove that in every iteration we
% can find a g(k)∈ (a(k),b(k)) where a(k) and b(k) are some constants not
% greater than 1 and that g(k) is the step which ensures that the recursive
% method of finding the minimum of the function f is acceptable. Since
% g(k)∈ (a(k),b(k)) , on inner optimization problem we look for g(k) in
% interval (0,1). I will use Golden Section Method and that
% method will find an interval for g(k), so in every (kth) iteration I will
% have an estimation for the optimal g(k) I search.

% This function needs, as inputs, the function f we want to minimize in the
% main problem, the current searching point (i.e. the integer k which is a
% counter of iterations of the main algorithm and the vector x,y from the
% main algorithm) and the vector d which is decided by the main method we
% use each time. Returns, as output, the estimation for descent step g that is the solution of inner optimization. 


function DescentStep = DescentStepByInnerOptimization(f,x,y,k,d)

% Firstly, we define the function of inner ptimization problem which is
% function of descent step g.
% We also define the interval we are going to look for g, so this is the
% inteval (0,1) as we analyse before. Then we define the accuracy "l" of
% estimation of the final interval the step g will be in. We use a very small "l" in order to have a good estimation. 
% At last, we will choose the middle of that final interval as final
% estimation for g. 

Fg=@(g) f( x(k)+ g.*d(1), y(k) +g.*d(2) );

a=0;
b=1;
l=0.00001;

[~, ~ , FinalInterval, ~, ~] = GoldenSectionMethod(Fg,a,b,l);

DescentStep = (FinalInterval(1)+FinalInterval(2) ) / 2;


end