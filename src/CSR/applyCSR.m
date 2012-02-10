function [beam W dE z]=applyCSR(beam,ind,nbin,smoothVal)
% [beam W dE z]=applyCSR(beam,ind,nbin,smoothVal)
%  Calculate CSR wake for this element and apply to beam
%
% beam: Lucretia beam
% ind: BEAMLINE index
% nbin: number of bins to use for CSR calculation
% smoothVal: level of smoothing to use (int >= 1)
global BEAMLINE
W=[]; dE=[]; z=[];

% If zero length element just return
if ~isfield(BEAMLINE{ind},'L') || BEAMLINE{ind}.L==0
  return
end
% Flip z definition
beam.Bunch.x(5,:)=-beam.Bunch.x(5,:);
% Bin beam particle longitudinal direction
zmin=min(beam.Bunch.x(5,:));
zmax=max(beam.Bunch.x(5,:));
if zmin==zmax
  error('Need some spread in z-distribution of bunch to compute CSR!')
end
z=linspace(zmin,zmax,nbin);
[count,bininds] = histc(beam.Bunch.x(5,:),z);
q = accumarray(bininds',beam.Bunch.Q')';
z=z-mean(z);
bw=abs(z(2)-z(1));
Q=sum(q);
q=q./bw; %(q./Q)./bw;
q=smoothn(q,smoothVal);
dq=[diff(q)./bw 0];

% Electron charge
qe=-1.6021773e-19;

% Is this a bend or drift?
if strcmp(BEAMLINE{ind}.Class,'SBEN')
  % Find distance from start of bend
  if ind==1
    ind1=ind;
  else
    for iele=ind-1:-1:0
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
  R=BEAMLINE{ind}.L/(2*sin(abs(BEAMLINE{ind}.Angle)/2));
  PHI=0;
  for iele=ind1:ind
    if isfield(BEAMLINE{iele},'Angle')
      PHI=PHI+abs(BEAMLINE{iele}.Angle);
    end
  end
  SL=(R*PHI^3)/24;
  % Loop over particle distribution and form wakefield function and
  % calculate energy loss for each bin
  W=arrayfun(@(x) csrInt1(x,z,SL,bw,dq,R,PHI,q),1:nbin);
  dE=(W.*Q.*BEAMLINE{ind}.L)./(1e9*qe); % GeV
else % DRIFT or other element following bend
  % Get parameters for wake calculation
  % -- find bend elements upstream
  bendele=[];
  for iele=ind-1:-1:1
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
  phi=abs(sum(arrayfun(@(x) BEAMLINE{x}.Angle,bendele)));
  L=sum(arrayfun(@(x) BEAMLINE{x}.L,bendele));
  R=L/(2*sin(phi/2));
  lDecay=3*(24*std(beam.Bunch.x(5,:))*R^2)^(1/3);
  % --- distance from bend
  X=((BEAMLINE{ind}.S+BEAMLINE{ind}.L)-(BEAMLINE{bendele(1)}.S+BEAMLINE{bendele(1)}.L))/R;
  if X*R > lDecay
    beam.Bunch.x(5,:)=-beam.Bunch.x(5,:);
    return
  end
  dsmax=((R*phi^3)/24)*((phi+4*X)/(phi+X));
  % Pre-compute psi function
  psi=zeros(1,length(z)-1);
  for isp=2:length(z)
    ds=abs(z(1)-z(isp));
    a=24*ds/R; b=4*X; C=[-1 -b 0 a a*X];
    rpsi=roots(C);
    psi(isp-1)=max(rpsi(arrayfun(@(x) ~any(imag(x)),rpsi)));
  end
  % Calculate CSR wake and energy loss per bin
  W=arrayfun(@(x) csrInt2(x,z,dsmax,psi,X,dq,q,phi,bw,R),1:nbin);
  dE=(W.*Q.*BEAMLINE{ind}.L)./(1e9*qe); % GeV
end
% Apply energy loss for all particles in each bin
beam.Bunch.x(6,:)=beam.Bunch.x(6,:)+dE(bininds);
beam.Bunch.x(5,:)=-beam.Bunch.x(5,:);

function W=csrInt1(is,z,SL,bw,dq,R,PHI,q)
isp=find(z>(z(is)-SL),1,'first'):is-1;
ZINT=sum((1./(z(is)-z(isp)).^(1/3)).*dq(isp).*bw);
[Y I1]=min(abs(z-(z(is)-((R*PHI^3)/6))));
[Y I2]=min(abs(z-(z(is)-((R*PHI^3)/24))));
W=-(4/(R*PHI))*q(I1) + (4/(R*PHI))*q(I2) + (2/((3*R^2)^(1/3)))*ZINT;

function W=csrInt2(is,z,dsmax,psi,X,dq,q,phi,bw,R)
isp=find(z>(z(is)-dsmax),1,'first'):is-1;
ZINT=sum((1./(psi(is-isp)+2.*X)).*dq(isp).*bw);
[Y I]=min(abs(z-(z(is)-(R/6)*phi^2*(phi+3*X))));
[Y I1]=min(abs(z-(z(is)-dsmax)));
W = (4/R)*( (q(I1)/(phi+2*X)) + ZINT ) - (4/R)*(1/(phi+2*X))*q(I) ;