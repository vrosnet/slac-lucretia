function f=gengauss_min(p, arg1)

%  f=gengauss_min(p, arg1);
%
%  Returns the error between the data and the values computed by the current
%  function of p.  Assumes a function of the form:
%
%  y=c(1)+c(2)*(1-exp(-abs((x-p(1))/(p(2)*sqrt(gamma(1/p(3))/gamma(3/p(3)))))^p(3)))
%
%  with 2 linear parameters and 3 nonlinear parameters (see also gengauss).
%
%  INPUTS:   p:      A vector of 3 scalar asymmetric gaussian parameters:
%
%                    p(1) => Horizontal offset (<x>).
%                    p(2) => Standard deviation (sigma) of the generialized 
%                            gaussian.
%                    p(3) => Shape parameter (=2 for symmetric gaussian).
%
%            arg1:   x (arg1(:,1)), y (arg1(:,2)), and dy (arg1(:,3)) to fit.
%
%  OUTPUTS:  f:      Computed error value.

%===============================================================================

x=arg1(:,1);
y=arg1(:,2);
dy=arg1(:,3);

A=[ones(size(x)) gengauss(x,p(1),p(2),p(3))];
c=A\y;
z=A*c;
f=norm((z-y)./dy);

