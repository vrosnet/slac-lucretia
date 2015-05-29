function [yfit,q,dq,chisq]=gengauss_fit(x,y,dy,y_off)

%  [yfit,q,dq,chisq]=gengauss_fit(x,y[,dy,y_off]);
%
%  Non-linear fitting routine to fit an  generalized gaussian "bell" curve
%  to the data in vectors "x" and "y", where "x" contains the independent
%  variable data and "y" contains the dependent data.  The fit form is as
%  follows:
%
%  yfit=q(1)+q(2)*exp(-abs((x-q(3))/(q(4)*sqrt(gamma(1/q(5))/gamma(3/q(5)))))^q(5))
%
%  When "y_off"=0 then q(1) is returned as exactly zero.
%
%  INPUTS:   x:	     A vector (row or column) of independent variable
%	     	     data.
%	     y:	     A vector (row or column) of dependent data.
%	     dy:     (Optional,DEF=no errors) A vector (row or column) of
%                    the errors on the dependent variable data.
%            y_off:  (Optional,DEF=1) If "y_off" .NE. 0, then the gaussian
%                    fit includes a "y" offset parameter in the fit.
%                    If "y_off" = 0, this forces the fit to include no
%                    "y" offset and therefore q(1) is set exactly to zero.
%
%  OUTPUTS:  yfit:   A vector of fitted bell curve data from which the
%                    difference to the original data "y" has been minimized.
%            q:      A vector of the 5 scalars which are fitted:
%
%                    q(1) => DC offset in the data (if "y_off"=0 then q(1)=0).
%                    q(2) => Amplitude scaling factor of asymmetric gaussian.
%                    q(3) => Horizontal offset (<x>).
%                    q(4) => Distribution width of the generalized gaussian.
%			     Corresponds to gaussian standard deviation.
%			     q(4)^2 = variance.	
%                    q(5) => "bell" shape parameter (=2 for gaussian).
%
%            dq:     A vector of the 5 errors on each of the above q's.
%            chisq:  Chisquare per degree of freedom (~1 => good fit when
%                   "dy" is given)

%===============================================================================
 
if length(y)<6
  error('Need at least 6 data points to fit an asymmetric gaussian')
end

x=x(:);
y=y(:);
if ~exist('dy','var'),dy=ones(size(y));end
if isempty(dy),dy=ones(size(y));end
dy=dy(:);
if ~exist('y_off','var'),y_off=1;end

arg1(:,1)=x;
arg1(:,2)=y;
arg1(:,3)=dy;

% compute initial guesses for <x>, sigma, and shape parameter

[~,imax]=max(abs(y));
x_ymax=x(imax);
p=[x_ymax std(x) 2];

%p=fminsearch('agauss_min',p,[],arg1);
p=fminsearch('gengauss_min',p,[],arg1);
%z=sqrt(2*pi)*p(2)*agauss(x,p(1),p(2),p(3));
z=(2*sqrt(gamma(1/p(3))/gamma(3/p(3)))*p(2)*gamma(1+1/p(3)))*gengauss(x,p(1),p(2),p(3));
if y_off
  Q=[ones(size(z)) z];
  [yfit,~,c]=fit(Q,y,dy);
  q=[c(1) c(2) p(1) p(2) p(3)];
else
  Q=z;
  [yfit,~,c]=fit(Q,y,dy);
  q=[0 c(1) p(1) p(2) p(3)];
end
chisq=norm((y-yfit)./dy)/sqrt(length(y)-length(q));
chisq=chisq^2;
%covm=agauss_cov(x,dy,q(1),q(2),q(3),q(4),q(5));
%da=sqrt(abs(covm(1,1)));
%db=sqrt(abs(covm(2,2)));
%dc=sqrt(abs(covm(3,3)));
%dd=sqrt(abs(covm(4,4)));
%de=sqrt(abs(covm(5,5)));
%dq=[da db dc dd de];
%q(4)=abs(q(4));

%Set to 0. Will be changed if everything works well.
da=0;
db=0;
dc=0;
dd=0;
de=0;
dq=[da db dc dd de];
q(4)=abs(q(4));
