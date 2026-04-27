% The Armijo Step Size Rule is a successive reduction method of the step
% g(k) which aims to satisfy the criterion: 
% f(x(k+1),y(k+1))<=f(x(k),y(k))+a.*g(k).*( dot( transpose(d), Gradientf(x(k),y(k)) ) )

% This function needs (as inputs) the function f, the current searching point
% [x(k),y(k)], the parameters a and b, the initial step s and the vector d of kth iteration of the main algorithm.
% In this method we choose the step as g(k)=s.*(b.^m(k)) and we try to find
% the minimum non-negative integer m(k) (of the kth iteration) which satisfies the criterion above.
% Returns, as outputs, that integer m(k) and the resulting step g(k)= s.*(b.^m(k)).

% THe Armijo Step Size Rule is a method which is used to determine the step
% g(k) we need to use in an another algorithm.
% !!! Sign: the integer k regards to the iterations of the main algorithm, not the iterations of the Armijo Method.
% So k is a constant for the Armijo Method. 

% !! The parameters a and b: They usually take values like: a∈[0.00001 , 0.1]
% and b∈[0.1 , 0.5].

function [DescentStep, IntegerM] = ArmijoStepSizeRule(f,s,a,b,d,x,y,k)
Gradientf=gradient(f);

% The integer m(k), which we are looking for, is the minimum non-negative integer
% which satisfies the criterion:
% f(x(k),y(k))-f(x(k+1),y(k+1))>= -a.*g(k).*( transpose(d).*Gradientf(x(k),y(k)) )
% where g(k)= s.*(b.^m(k)). 
% So, firstly we set m(k)=0 and in every iteration we need to determine the
% x(k+1) and y(k+1) using the respective m(k) (i.e. the respective g(k)).

m(k)=0;
g(k)=s(1).*(b(1).^m(k));
x(k+1)=x(k)+ g(k).*d(1);
y(k+1)=y(k)+ g(k).*d(2);

% We set two quantities A(k) and B(k) for simplicity.
A(k)=f(x(k),y(k))-f(x(k+1),y(k+1));
B(k)=-a(1).*(b(1).^m(k)).*( dot( transpose(d), Gradientf(x(k),y(k)) ) );

while(A(k)<B(k))
m(k)=m(k)+1;

g(k)=s.*(b.^m(k));
x(k+1)=x(k)+ g(k).*d(1);
y(k+1)=y(k)+ g(k).*d(2);

% Inform the quantities A(k) and B(k).
A(k)=f(x(k),y(k))-f(x(k+1),y(k+1));
B(k)=-a(1).*(b(1).^m(k)).*( dot( transpose(d), Gradientf(x(k),y(k)) ) );
end

% Return the values
DescentStep=g(k);
IntegerM=m(k);

end