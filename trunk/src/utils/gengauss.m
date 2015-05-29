function y=gengauss(x,x_bar,sig,p)

%  y=gengauss(x[,x_bar,sig,p]);
%
%  Generalized gaussian function to create a "bell" curve which shape is
%  controlled by the parameter p. In particular case when p=2 the function 
%  becomes a well known gaussian function, becomes more narrow when p<2 and more
%  broad when p>2. The function is normalized so that its variance is equal to 
%  sigma^2 and the area under the full curve will be 1.
%  gamma(a,x) denotes complementary incomplete gamma function defined by
%  gamma(a,x) = int_{x}^{infty} t^{a-1} e^{-t} dt. 
%
%
%  y=exp(-abs((x-x_bar)/(sig*A(p)))^p) / (2*A(p)*sig*gamma(1+1/p))
%  A(p)=sqrt(gamma(1/p)/gamma(3/p))
% 
%  for p=2: A(2)=sqrt(2), 2*gamma(3/2) = sqrt(2) and
%  y=exp(-abs((x-x_bar)/(sig*sqrt(2)))^2) / (sqrt(2*pi)*sig) -> Gaussian funct.
%
%  INPUTS:   x:      The independent variable in the exponential of the
%                    gaussian (column or row vector)
%            x_bar:  (Optional,DEF=0) The center of the curve on the
%                    "x" axis (mean) which defaults to 0 if not given
%                    (scalar)
%            sig:    (Optional,DEF=1) The distribution width, sig^2 = variance.
%                     Default is 1 if sig is not given (scalar)
%            p:      (Optional,DEF=2) The shape parameter
%                    which defaults to 2 if not given (scalar)
%
%  OUTPUTS:  y:      The values of the generalized gaussian at each "x"
%                    (vector the same size as "x").
%
%  EXAMPLE:  >> x=[-2.5:.5:2.5];
%            >> y=gengauss(x);
%            >> y
%							y =
%						  Columns 1 through 7
%					    0.0175    0.0540    0.1295    0.2420    0.3521    0.3989    0.3521
%						  Columns 8 through 11
%    					0.2420    0.1295    0.0540    0.0175

%===============================================================================
if nargin==1
  x_bar=0;
  sig=1;
  p=2;
elseif nargin==2
  sig=1;
  p=2;
elseif nargin==3
  p=2;
end

if sig==0
   sigx=1e-3;
else
   sigx=sig;
end

y=exp(-abs((x-x_bar)/(sqrt(gamma(1/p)/gamma(3/p))*sigx)).^p)/(2*sqrt(gamma(1/p)/gamma(3/p))*sigx*gamma(1+1/p));

