        function y = gaussHaloFn(x,x_bar,sig,a,b,haloSIGmin)
              
%       y = gaussHaloFn(x,x_bar,sig,a,b,haloSIGmin);
%
%       Gaussian function to create a "bell" curve of width "sig" centered
%       around "x_bar" from the independent variable "x".  It is normalized
%       already so that the area under the full curve will be 1.
%       + Halo distribution parameterization as a*x^-b (x>haloXmin)
%
%                                             2      2       -b
%               y = (1/sqrt(2pi)) exp(-(x-<x>) /2 sig ) + a.x
%
%     INPUTS:   x:      The independent variable in the exponential of the
%                       gaussian (column or row vector)
%               sig:    The gaussian standard deviation 
%                       ("width") which defaults to 1 if not given (scalar)
%               x_bar:  The center of the gaussian on the 
%                       "x" axis (mean) which defaults to 0 if not given
%                       (scalar)
%               a,b,haloSIGmin: power-law fit to halo a.x^-b (> haloSIGmin)
%
%     OUTPUTS:  y:      The values of the fit at each "x" (vector the
%                       same size as "x").

%===============================================================================

if ~exist('sig','var')
  sig = 1;
end
if ~exist('x_bar','var')
  x_bar = 0;
end
if ~exist('a','var')
  a=0;
end
if ~exist('b','var')
  b=1;
end
if ~exist('haloSIGmin','var')
  haloSIGmin=1;
end
if sig<1e-3
  sig=1e-3;
end
y = (1/(sqrt(2*pi)*sig))*exp(-( (x-x_bar).^2)/(2*sig^2)) + ...
  double((abs(x)./sig)>haloSIGmin).*a.*abs(x-x_bar).^-b;
