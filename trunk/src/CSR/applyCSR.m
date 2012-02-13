function [beam W dE zOut]=applyCSR(beam,nbin,smoothVal,itrack)
% [beamP W dE z]=applyCSR(beam,beamQ,ind,nbin,smoothVal,X)
%  Calculate CSR wake and change provided momentum profile
%
% beamZ: Lucretia Bunch.x
% beamQ: charge of macro particles
% nbin: number of bins to use for CSR calculation
% smoothVal: level of smoothing to use (int >= 1)
global BEAMLINE
persistent lastsz bininds z Z ZSP
W=[]; dE=[]; zOut=[];

%- If zero length element just return
if ~isfield(BEAMLINE{itrack},'L') || BEAMLINE{itrack}.L==0
  return
end

% Find distance from start of bend
if strcmp(BEAMLINE{itrack}.Class,'SBEN')
  if itrack==1
    ind1=itrack;
  else
    for iele=itrack-1:-1:0
      % is this still a bend?
      if iele==0 || ~strcmp(BEAMLINE{iele}.Class,'SBEN')
        % is it a real element that isn't a bend?
        % - then the start of bend is the element after this
        if iele==0 || (isfield(BEAMLINE{iele},'L') && BEAMLINE{iele}.L>0)
          ind1=iele+1;
          break
        end
      end
    end
  end
  % Get parameters for CSR wake calculation
  R=BEAMLINE{itrack}.L/(2*sin(abs(BEAMLINE{itrack}.Angle)/2));
  PHI=0;
  for iele=ind1:itrack
    if isfield(BEAMLINE{iele},'Angle')
      PHI=PHI+abs(BEAMLINE{iele}.Angle);
    end
  end
  X=0;
else
  % -- find bend elements upstream
  bendele=[];
  for iele=itrack-1:-1:1
    if strcmp(BEAMLINE{iele}.Class,'SBEN')
      bendele(end+1)=iele;
    elseif ~isempty(bendele) && isfield(BEAMLINE{iele},'L') && BEAMLINE{iele}.L>0
      break
    end
  end
  if isempty(bendele)
    error('No BEND found upstream for CSR calculation')
  end
  % --- bend angle and radius
  PHI=abs(sum(arrayfun(@(x) BEAMLINE{x}.Angle,bendele)));
  L=sum(arrayfun(@(x) BEAMLINE{x}.L,bendele));
  R=L/(2*sin(PHI/2));
  % --- distance from bend
  X=((BEAMLINE{itrack}.S+BEAMLINE{itrack}.L)-(BEAMLINE{bendele(1)}.S+BEAMLINE{bendele(1)}.L))/R;
  lDecay=3*(24*std(beam.x(5,:))*R^2)^(1/3);
  if X*R > lDecay; return; end;
end

% Generate longitudinal grid, only when beam length changes by
% more than 10%
if isempty(z) || (abs(std(beam.x(5,:))-lastsz)/lastsz)>0.1
  zmin=min(-beam.x(5,:));
  zmax=max(-beam.x(5,:));
  if zmin==zmax
    error('Need some spread in z-distribution of bunch to compute CSR!')
  end
  z=linspace(zmin,zmax,nbin);
  [count,bininds] = histc(-beam.x(5,:),z);
  [Z ZSP]=meshgrid(z,z);
  lastsz=std(beam.x(5,:));
end

% Bin beam particle longitudinal direction
q = accumarray(bininds',beamQ')';
bw=abs(z(2)-z(1));
Q=sum(q);
q=q./bw; %(q./Q)./bw;
q=smoothn(q,smoothVal);
dq=[diff(q)./bw 0];

% Electron charge
qe=-1.6021773e-19;

% Is this a bend or drift?
if strcmp(BEAMLINE{itrack}.Class,'SBEN')
  SL=(R*PHI^3)/24;
  % Loop over particle distribution and form wakefield function and
  % calculate energy loss for each bin
  ZINT=zeros(nbin,1);
  for is=1:nbin
    isp=z>(z(is)-SL) & (1:length(z))<is ;
    ZINT(is)=sum((1./(z(is)-z(isp)).^(1/3)).*dq(isp).*bw);
  end
  IND1=abs(Z-(ZSP-(R*PHI^3)/6));
  IND2=abs(Z-(ZSP-(R*PHI^3)/24));
  [Y I1]=min(IND1,[],2);
  [Y I2]=min(IND2,[],2);
  W=-(4/(R*PHI)).*q(I1)' + (4/(R*PHI)).*q(I2)' + (2/((3*R^2)^(1/3))).*ZINT;
  dE=(W'.*Q.*BEAMLINE{itrack}.L)./(1e9*qe); % GeV
else % DRIFT or other element following bend
  % Get parameters for wake calculation
  dsmax=((R*PHI^3)/24)*((PHI+4*X)/(PHI+X));
  % Pre-compute psi function
  psi=zeros(1,length(z)-1);
  for isp=2:length(z)
    ds=abs(z(1)-z(isp));
    a=24*ds/R; b=4*X; C=[-1 -b 0 a a*X];
    % Polynomial roots via a companion matrix
    a = diag(ones(1,3),-1);
    d = C(2:end)./C(1);
    a(1,:) = -d;
    rpsi = eig(a);
    psi(isp-1)=max(rpsi(imag(rpsi)==0)); % take real root > 0
  end
  % Calculate CSR wake and energy loss per bin
  ZINT=zeros(nbin,1);
  for is=1:nbin
    isp=z>(z(is)-dsmax) & (1:length(z))<is ;
    ZINT(is)=sum((1./(psi(isp)+2.*X)).*dq(isp).*bw);
  end
  IND1=abs(Z-(ZSP-(R/6)*PHI^2*(PHI+3*X)));
  [Y I]=min(IND1,[],2);
  IND2=abs(Z-(ZSP-dsmax));
  [Y I1]=min(IND2,[],2);
  W = (4/R)*( (q(I1)'./(PHI+2*X)) + ZINT ) - (4/R)*(1/(PHI+2*X)).*q(I)' ;
  dE=(W'.*Q.*BEAMLINE{itrack}.L)./(1e9*qe); % GeV
end
% Apply energy loss for all particles in each bin
beam.x(6,:)=beam.x(6,:)+dE(bininds);

zOut=z;
