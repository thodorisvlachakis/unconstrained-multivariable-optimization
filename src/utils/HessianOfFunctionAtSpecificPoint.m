% This function returns the hessian of a function (with two variables)
% calculated at a given point.

function Hessian= HessianOfFunctionAtSpecificPoint(f,x0)
syms x y;

A(x,y)=hessian(f);

Hessian=A(x0(1),x0(2));


end