function [fitTerm,fitCoef,bsize_corrected,bsize,p] = beamTerms(dim,beam,doGaussfit)
% [fitTerm,fitCoef,bsize_corrected,bsize] = beamTerms(dim,beam [,doGaussfit])
% Show relative contributions of up to 3rd-order beam correlations to the
% beam size in the requested dimension
% 3rd order polynomial fit to 4 independent variables vs specified (dim)
% dependent variable
% REQUIRES: PolyfitnTools directory in search path
% INPUTS:
% -------
% dim = required lucretia beam dimension to correlate
% beam = lucretia beam
% doGaussfit (optional): if provided and true then bsize returned parameter
% shows effect on gaussian fit to core of distribution, else it is just the
% RMS
%
% OUTPUTS:
% --------
% fitTerm(35,6) = Beam correlation (e.g. if dim=3 [1 0 0 0 0 0] = SIGMA(3,1)
%                                                 [1 0 0 1 0 0] = T314
%                                                 [2 1 0 0 0 0] = U3112 etc... )
% fitCoef(1,35) = coefficient related to corresponding beam correlation
% bsize_corrected(1,3) = RMS beam size in dimension dim when all 1st-3rd order
%                        correlations are subtracted
% bsize(1,35) = RMS beam size improvement due to corresponding correlation term
%               being subtracted from beam
% p = results structure from polyfitn

% Warnings not caring about here
warnstate=warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

if ~exist('doGaussfit','var'); doGaussfit=false; end;

% remove constant offsets
for idim=1:6
  beam.Bunch.x(idim,:)=beam.Bunch.x(idim,:)-mean(beam.Bunch.x(idim,:));
end

% Perform 1st-3rd order fit
allpar=[1 2 3 4 5 6];
xfit=beam.Bunch.x(allpar(~ismember(allpar,dim)),:)';
yfit=beam.Bunch.x(dim,:)';
if doGaussfit
  nbin=max([length(yfit)/100 100]);
  [ fx , bc ] = hist(yfit,nbin) ;
  [~, q] = gauss_fit(bc,fx) ;
  bsize_initial = q(4) ;
else
  bsize_initial = std(yfit) ;
end
bsize_corrected=zeros(1,3);
for iorder=1:3
  p=polyfitn(xfit,yfit,iorder);
  bsize_corrected(iorder)=std(yfit-polyvaln(p,xfit));
end
xfit=beam.Bunch.x';
p=polyfitn(xfit,yfit,iorder);
% Get contribution of beamsize to each fit term
[Y I]=sort(sum(p.ModelTerms,2));
bsize=zeros(length(p.Coefficients),1);
fitTerm=zeros(length(p.Coefficients),length(p.ModelTerms(1,:)));
fitCoef=zeros(1,length(p.Coefficients));
for iterm=1:length(p.Coefficients)
  p2=polyfitn(xfit,yfit,p.ModelTerms(I(iterm),:));
  if Y(iterm)>0 && ~(length(find(p.ModelTerms(I(iterm),:)))==1 && find(p.ModelTerms(I(iterm),:))==dim)
    yfit=yfit-polyvaln(p2,xfit);
    if doGaussfit
      nbin=max([length(yfit)/100 100]);
      [ fx , bc ] = hist(yfit,nbin);
      [~, q] = gauss_fit(bc,fx) ;
      bsize(iterm)=abs(bsize_initial-q(4));
      bsize_initial=q(4);
    else
      bsize(iterm)=abs(bsize_initial-std(yfit));
      bsize_initial=std(yfit);
    end
  end
  fitTerm(iterm,:)=p.ModelTerms(I(iterm),:);
  fitCoef(iterm)=p2.Coefficients;
end

% Put warnings back to original state
warning(warnstate.state,'MATLAB:nearlySingularMatrix');