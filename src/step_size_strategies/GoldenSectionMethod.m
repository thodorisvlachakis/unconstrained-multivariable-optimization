% In this method we want the function f(x) and the interval [a,b} to be
% defined. Also, we want the same for the accuracy l. So we give them as
% inputs in our function.

function [A, B , FinalInterval, Iterations, Calculations]= GoldenSectionMethod(f,a,b,l)
% In the first iteration we choose the interval [a,b]. Also, we define the
% constant 0.618 which is the golden section. In every iteration we choose
% the searching points x1(k) x2(k).
% the variable k will be a counter for the iterations of the algorithm
% the value calc will be a counter for the calculations of function f(x).

g=0.618;
k=1;
calc=0;


a(k)=a;
b(k)=b;

x1(k)=a(k) + (1-g)*(b(k)-a(k));
x2(k)=a(k) + g*(b(k)-a(k));



while((b(k)-a(k))>l)
    
    if(f(x1(k))>f(x2(k)))
        % The new interval will be the (x1k,bk] and we also choose the next
        % searching points as x1(k+1)=x2(k) and x2(k+1)=a(k+1) + g*(b(k+1)-a(k+1)).

        a(k+1)=x1(k);
        b(k+1)=b(k);
        x1(k+1)=x2(k);
        x2(k+1)=a(k+1) + g*(b(k+1)-a(k+1));
    else
        % The new interval will be the [ak,x2k] and we also choose the next
        % searching points as x2(k+1)=x1(k) and x1(k+1)=a(k+1) + (1-g)*(b(k+1)-a(k+1)).
        
        a(k+1)=a(k);
        b(k+1)=x2(k);
        x2(k+1)=x1(k);
        x1(k+1)=a(k+1) + (1-g)*(b(k+1)-a(k+1));
        
    end

    % We inform thw counter of iterations of alorithm and the counter of
    % calculations of function f.in this method we already know one of
    % the values f(x1(k+1)) or f(x2(k+1)) because in every iteration we set the new searching points as x1(k+1)=x2(k) or x2(k+1)=x1(k).
    % That's why the counter of calculations increases by one on each
    % iteration. However, if k=1 we need two calculation because this is
    % the first iteration.
    if(k==1)
        calc=1;
    end

    k=k+1;
    calc=calc+1;

end

% I save the subintervals on each iteration in two vectors a and b and I
% give them back as an output.
A=a;
B=b;
FinalInterval=[a(k) b(k)];

Iterations=k;
Calculations=calc;


end