function [stat,emitData] = emit2dOTR(otruse,dointrinsic,printData)
% [stat emitData] = emit2dOTR(otruse,[dointrinsic,printData])
% Compute beam sigma matrix, covariance matrix, and chisquare from measured
% beam sizes and their errors at 3 or more OTRs
% (assumes data already gathered by OTR system and data input into EPICS PVs)
%
% Inputs:
% otruse = [1 * 4] logical array of which OTRs to use
% dointrinsic = logical (optional)
%               true = show intrinsic emittance data
%               false = show projected emittance data [default]
% printData = logical (optional)
%             true = display emittance data to screen and show plots
%             false = no printing to screen or plotting [default]
% Outputs:
% stat = Lucretia status return
% emitData = data formatted for sending over FlECS interface:
%   [energy ...
%    emitx demitx emitxn demitxn embmxn dembmxn ...
%    bmagx dbmagx bcosx dbcosx bsinx dbsinx ...
%    betax dbetax bx0 alphx dalphx ax0 chi2x ...
%    emity demity emityn demityn embmyn dembmyn ...
%    bmagy dbmagy bcosy dbcosy bsiny dbsiny ...
%    betay dbetay by0 alphy dalphy ay0 chi2y ...
%    ido length(id) id length(S) S sigxf sigx dsigx sigyf sigy dsigy ...
%    R{1:notr} ...
%    exip0 bxip0 axip0 eyip0 byip0 ayip0 ...
%    sigxip dsigxip sigpxip dsigpxip betaxip dbetaxip alphxip dalphxip ...
%    sigyip dsigyip sigpyip dsigpyip betayip dbetayip alphyip dalphyip ...
%   ];
% ------------------------------------------------------------------------------
% 04-Nov-2012, M. Woodley
%    Use energy=FL.SimModel.Initial.Momentum (BH1R fudge updated 24Oct12);
%    add to emitData: 4x4 OTR-to-OTR Rmats, design IP Twiss propagated IP beam
%    parameters
% ------------------------------------------------------------------------------

global BEAMLINE INSTR FL PS %#ok<NUSED>
stat{1}=1;
emitData=[];

debug=0;

% Parse optional parameters
if ~exist('dointrinsic','var') || isempty(dointrinsic)
  dointrinsic=false;
end
if ~exist('printData','var') || isempty(printData)
  printData=false;
end

% which OTRs to use
if sum(otruse)<3
  stat{1}=-1; stat{2}='Must select 3 or 4 OTRs for use';
  return
end
oname={};
for iotr=find(otruse)
  oname{end+1}=sprintf('OTR%dX',iotr-1);
end
notr=length(oname);
z=zeros(size(oname));
sigx=z;dsigx=z;sigy=z;dsigy=z;theta=z;

% get pointers
ido=z;idoi=z;
for n=1:notr
  ido(n)=findcells(BEAMLINE,'Name',oname{n});
  idoi(n)=findcells(INSTR,'Index',ido(n));
end

% Get OTR data from PVs
if (debug)
  load userData/emit2dOTRtest %#ok<*UNRCH>
  id=find(otruse);
  for n=1:notr
    if dointrinsic
      sigx(n)=1e-6*sqrt(pvdata(1,id(n))); % m
      dsigx(n)=1e-6*sqrt(pvdata(2,id(n))); % m
      sigy(n)=1e-6*sqrt(pvdata(3,id(n))); % m
      dsigy(n)=1e-6*sqrt(pvdata(4,id(n))); % m
    else
      sigx(n)=1e-6*pvdata(7,id(n)); % m
      dsigx(n)=1e-6*pvdata(8,id(n)); % m
      sigy(n)=1e-6*pvdata(9,id(n)); % m
      dsigy(n)=1e-6*pvdata(10,id(n)); % m
    end
  end
  z=zeros(1,notr);
  theta=z; % deg
  DX=z;DPX=z;DY=z;DPY=z; % m
  err=1e-6*ones(1,notr);
  dDX=err;dDPX=err;dDY=err;dDPY=err; % m
else
  for n=1:notr
    otrNum=str2double(oname{n}(4))+1;
    [stat,data]=getOTRsize(otrNum);
    ictdata=data.ict; icterrdata=data.icterr; %#ok<NASGU>
    if (stat{1}~=1)
      stat{1}=-1;
      stat{2}=sprintf('Failed to get data for OTR%dX',n-1);
      return
    end
    rawotrdata{n}=data; %#ok<NASGU>
    if dointrinsic
      sigx(n)=1e-6*sqrt(data.sig11);dsigx(n)=1e-6*(data.sig11err/2); % m
      sigy(n)=1e-6*sqrt(data.sig33);dsigy(n)=1e-6*(data.sig33err/2); % m
    else
      sigx(n)=1e-6*data.projx;dsigx(n)=1e-6*data.projxerr; % m
      sigy(n)=1e-6*data.projy;dsigy(n)=1e-6*data.projyerr; % m
      %sigx(n)=1e-6*data.sigx;dsigx(n)=1e-6*data.projxerr; % m
      %sigy(n)=1e-6*data.sigy;dsigy(n)=1e-6*data.projyerr; % m
    end
    theta(n)=data.theta; % deg
  end
  % get dispersion data
  DX=z;DPX=z;DY=z;DPY=z;dDX=z;dDPX=z;dDY=z;dDPY=z;
  for n=1:notr
    D=INSTR{idoi(n)}.dispref;dD=INSTR{idoi(n)}.dispreferr;
    [DX(n),DPX(n),DY(n),DPY(n)]=deal(D(1),D(2),D(3),D(4)); % m,rad,m,rad
    [dDX(n),dDPX(n),dDY(n),dDPY(n)]=deal(dD(1),dD(2),dD(3),dD(4)); % m,rad,m,rad
  end
end
if isfield(FL,'props') && isfield(FL.props,'dE')
  dp=FL.props.dE;
else
  dp=8e-4; % nominal energy spread
end

% correct measured spot sizes for dispersion
sigxt=sigx;sigxd=abs(dp*DX);sigx2=sigxt.^2-sigxd.^2;
sigyt=sigy;sigyd=abs(dp*DY);sigy2=sigyt.^2-sigyd.^2;
if (any(sigx2<0)||any(sigy2<0))
  stat{1}=-1;
  stat{2}='Negative sigx2 or sigy2 values after dispersion correction';
  return
end
sigx=sqrt(sigx2); % dispersion corrected
sigy=sqrt(sigy2); % dispersion corrected
if dointrinsic
  txt{1}='Ellipse:';
else
  txt{1}='Projected:';
end
txt{2}='  sigxt   sigxd    sigx   sigyt   sigyd    sigy';
txt{3}='------- ------- ------- ------- ------- -------';
%       nnnn.nn nnnn.nn nnnn.nn nnnn.nn nnnn.nn nnnn.nn
for n=1:notr
  txt{3+n}=sprintf('%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f', ...
    1e6*[sigxt(n),sigxd(n),sigx(n),sigyt(n),sigyd(n),sigy(n)]);
end

% get the model
R=cell(1,notr);
Rx=zeros(notr,2);Ry=zeros(notr,2);
for n=1:notr
  [stat,Rab]=RmatAtoB(ido(1),ido(n));
  if (stat{1}~=1),error(stat{2}),end
  R{n}=Rab(1:4,1:4);
  Rx(n,:)=[Rab(1,1),Rab(1,2)];
  Ry(n,:)=[Rab(3,3),Rab(3,4)];
end

% get design Twiss at first OTR
energy=FL.SimModel.Initial.Momentum;
egamma=energy/0.51099906e-3;
ex0=FL.SimModel.Initial.x.NEmit/egamma;
bx0=FL.SimModel.Design.Twiss.betax(ido(1));
ax0=FL.SimModel.Design.Twiss.alphax(ido(1));
ey0=FL.SimModel.Initial.y.NEmit/egamma;
by0=FL.SimModel.Design.Twiss.betay(ido(1));
ay0=FL.SimModel.Design.Twiss.alphay(ido(1));

% get design Twiss at IP
ip=findcells(BEAMLINE,'Name','IP');
ex0ip=ex0;
bx0ip=FL.SimModel.Design.Twiss.betax(ip);
ax0ip=FL.SimModel.Design.Twiss.alphax(ip);
ey0ip=ey0;
by0ip=FL.SimModel.Design.Twiss.betay(ip);
ay0ip=FL.SimModel.Design.Twiss.alphay(ip);

% load analysis variables
dsigxd=dDX.*dp;
dsigyd=dDY.*dp;
x=sigx.^2;% m^2
dx=sqrt(((2.*sigxt).^2.*dsigx.^2)+((2.*sigxd).^2.*dsigxd.^2));
y=sigy.^2; % m^2
dy=sqrt(((2.*sigyt).^2.*dsigy.^2)+((2.*sigyd).^2.*dsigyd.^2));
x=x';dx=dx';y=y';dy=dy'; % columns

% compute least squares solution
Mx=zeros(notr,3);My=zeros(notr,3);
for n=1:notr
  Mx(n,1)=Rx(n,1)^2;
  Mx(n,2)=2*Rx(n,1)*Rx(n,2);
  Mx(n,3)=Rx(n,2)^2;
  My(n,1)=Ry(n,1)^2;
  My(n,2)=2*Ry(n,1)*Ry(n,2);
  My(n,3)=Ry(n,2)^2;
end
zx=x./dx;zy=y./dy;
Bx=zeros(notr,3);By=zeros(notr,3);
for n=1:notr
  Bx(n,:)=Mx(n,:)/dx(n);
  By(n,:)=My(n,:)/dy(n);
end
Tx=inv(Bx'*Bx);Ty=inv(By'*By);
u=Tx*Bx'*zx;chi2x=zx'*zx-zx'*Bx*Tx*Bx'*zx;du=sqrt(diag(Tx)); %#ok<NASGU,MINV>
v=Ty*By'*zy;chi2y=zy'*zy-zy'*By*Ty*By'*zy;dv=sqrt(diag(Ty)); %#ok<NASGU,MINV>

% convert fitted input sigma matrix elements to emittance, BMAG, ...
[p,dp]=emit_params(u(1),u(2),u(3),Tx,bx0,ax0);
p(1:3)=abs(p(1:3));
% if (any(imag(p(1:3))~=0)||any(p(1:3)<=0))
%   stat{1}=-1;
%   if dointrinsic
%     stat{2}='Error in horizontal intrinsic emittance computation';
%   else
%     stat{2}='Error in horizontal projected emittance computation';
%   end
%   return
% end
emitx=p(1);demitx=dp(1);
bmagx=p(2);dbmagx=dp(2);
embmx=p(3);dembmx=dp(3);
betax=p(4);dbetax=dp(4);
alphx=p(5);dalphx=dp(5);
bcosx=p(6);dbcosx=dp(6);
bsinx=p(7);dbsinx=dp(7);
emitxn=egamma*emitx;demitxn=egamma*demitx;
embmxn=egamma*embmx;dembmxn=egamma*dembmx;

[p,dp]=emit_params(v(1),v(2),v(3),Ty,by0,ay0);
p(1:3)=abs(p(1:3));
% if (any(imag(p(1:3))~=0)||any(p(1:3)<=0))
%   stat{1}=-1;
%   if dointrinsic
%     stat{2}='Error in vertical intrinsic emittance computation';
%   else
%     stat{2}='Error in vertical projected emittance computation';
%   end
%   return
% end
emity=p(1);demity=dp(1);
bmagy=p(2);dbmagy=dp(2);
embmy=p(3);dembmy=dp(3);
betay=p(4);dbetay=dp(4);
alphy=p(5);dalphy=dp(5);
bcosy=p(6);dbcosy=dp(6);
bsiny=p(7);dbsiny=dp(7);
emityn=egamma*emity;demityn=egamma*demity;
embmyn=egamma*embmy;dembmyn=egamma*dembmy;

% results txt
txt{end+1}=' ';
if dointrinsic
  txt{end+1}=sprintf('Horizontal intrinsic emittance parameters at %s',oname{1});
else
  txt{end+1}=sprintf('Horizontal projected emittance parameters at %s',oname{1});
end
txt{end+1}='-------------------------------------------------------';
txt{end+1}=sprintf('energy     = %10.4f              GeV',energy);
txt{end+1}=sprintf('emit       = %10.4f +- %9.4f pm',1e12*emitx,1e12*demitx);
txt{end+1}=sprintf('emitn      = %10.4f +- %9.4f nm',1e9*emitxn,1e9*demitxn);
txt{end+1}=sprintf('emitn*bmag = %10.4f +- %9.4f nm',1e9*embmxn,1e9*dembmxn);
txt{end+1}=sprintf('bmag       = %10.4f +- %9.4f      (%9.4f)',bmagx,dbmagx,1);
txt{end+1}=sprintf('bmag_cos   = %10.4f +- %9.4f      (%9.4f)',bcosx,dbcosx,0);
txt{end+1}=sprintf('bmag_sin   = %10.4f +- %9.4f      (%9.4f)',bsinx,dbsinx,0);
txt{end+1}=sprintf('beta       = %10.4f +- %9.4f m    (%9.4f)',betax,dbetax,bx0);
txt{end+1}=sprintf('alpha      = %10.4f +- %9.4f      (%9.4f)',alphx,dalphx,ax0);
txt{end+1}=sprintf('chisq/N    = %10.4f',chi2x);
txt{end+1}=' ';
if dointrinsic
  txt{end+1}=sprintf('Vertical intrinsic emittance parameters at %s',oname{1});
else
  txt{end+1}=sprintf('Vertical projected emittance parameters at %s',oname{1});
end
txt{end+1}='-------------------------------------------------------';
txt{end+1}=sprintf('energy     = %10.4f              GeV',energy);
txt{end+1}=sprintf('emit       = %10.4f +- %9.4f pm',1e12*emity,1e12*demity);
txt{end+1}=sprintf('emitn      = %10.4f +- %9.4f nm',1e9*emityn,1e9*demityn);
txt{end+1}=sprintf('emitn*bmag = %10.4f +- %9.4f nm',1e9*embmyn,1e9*dembmyn);
txt{end+1}=sprintf('bmag       = %10.4f +- %9.4f      (%9.4f)',bmagy,dbmagy,1);
txt{end+1}=sprintf('bmag_cos   = %10.4f +- %9.4f      (%9.4f)',bcosy,dbcosy,0);
txt{end+1}=sprintf('bmag_sin   = %10.4f +- %9.4f      (%9.4f)',bsiny,dbsiny,0);
txt{end+1}=sprintf('beta       = %10.4f +- %9.4f m    (%9.4f)',betay,dbetay,by0);
txt{end+1}=sprintf('alpha      = %10.4f +- %9.4f      (%9.4f)',alphy,dalphy,ay0);
txt{end+1}=sprintf('chisq/N    = %10.4f',chi2y);

% propagate measured beam to OTRs
xf=sqrt(Mx*u);yf=sqrt(My*v);
txt{end+1}=sprintf('Propagated spot sizes');
txt{end+1}='---------------------------------------------------------------------------';
for n=1:notr
  txt{end+1}=sprintf('%5s(x) = %6.1f um (%6.1f +- %6.1f), (y) = %6.1f um (%6.1f +- %6.1f)', ...
    oname{n},1e6*[xf(n),sigx(n),dsigx(n),yf(n),sigy(n),dsigy(n)]);
end
txt{end+1}=' ';

% back propagate fitted sigma matrices to MDISP
sigx1=[u(1),u(2);u(2),u(3)];sigy1=[v(1),v(2);v(2),v(3)];
id1=findcells(BEAMLINE,'Name','MDISP');
id2=findcells(BEAMLINE,'Name','BEGFF');
[stat,Rab]=RmatAtoB(id1,ido(1));
if (stat{1}~=1),error(stat{2}),end
Rx=Rab(1:2,1:2);Ry=Rab(3:4,3:4);
sigx0=inv(Rx)*sigx1*inv(Rx');sigy0=inv(Ry)*sigy1*inv(Ry'); %#ok<MINV>

% forward propagate through diagnostic section ...
S=FL.SimModel.Design.Twiss.S';
id=(id1:id2)';
idt=find([1;diff(S(id))]~=0); % unique S values
id=id(idt); %#ok<FNDSB>
sigxf=zeros(size(id));sigyf=zeros(size(id));
for n=1:length(id)
  [stat,Rab]=RmatAtoB(id1,id(n));
  if (stat{1}~=1),error(stat{2}),end
  Rx=Rab(1:2,1:2);Ry=Rab(3:4,3:4);
  sigxm=Rx*sigx0*Rx';sigym=Ry*sigy0*Ry';
  sigxf(n)=sqrt(sigxm(1,1));sigyf(n)=sqrt(sigym(1,1));
end

% propagate measured beam from first OTR to IP (ignore coupling)
[stat,Rab]=RmatAtoB(ido(1),ip); % first OTR to IP
if (stat{1}~=1),error(stat{2}),end

sig0=[u(1),u(2);u(2),u(3)];
Rx=Rab(1:2,1:2);
sigip=Rx*sig0*Rx';
RT=[   Rx(1,1)^2            2*Rx(1,1)*Rx(1,2)            Rx(1,2)^2   ; ...
    Rx(1,1)*Rx(2,1)  Rx(1,1)*Rx(2,2)+Rx(1,2)*Rx(2,1)  Rx(1,2)*Rx(2,2); ...
       Rx(2,1)^2            2*Rx(2,1)*Rx(2,2)            Rx(2,2)^2   ];
T=RT*Tx*RT'; % propagate covariance matrix to IP
[p,dp]=emit_params(sigip(1,1),sigip(1,2),sigip(2,2),T,1,0);
sigxip=sigip(1,1);dsigxip=sqrt(T(1,1));
sigpxip=sigip(2,2);dsigpxip=sqrt(T(3,3));
betaxip=p(4);dbetaxip=dp(4);
alphxip=p(5);dalphxip=dp(5);

sig0=[v(1),v(2);v(2),v(3)];
Ry=Rab(3:4,3:4);
sigip=Ry*sig0*Ry';
RT=[   Ry(1,1)^2            2*Ry(1,1)*Ry(1,2)            Ry(1,2)^2   ; ...
    Ry(1,1)*Ry(2,1)  Ry(1,1)*Ry(2,2)+Ry(1,2)*Ry(2,1)  Ry(1,2)*Ry(2,2); ...
       Ry(2,1)^2            2*Ry(2,1)*Ry(2,2)            Ry(2,2)^2   ];
T=RT*Ty*RT'; % propagate covariance matrix to IP
[p,dp]=emit_params(sigip(1,1),sigip(1,2),sigip(2,2),T,1,0);
sigyip=sigip(1,1);dsigyip=sqrt(T(1,1));
sigpyip=sigip(2,2);dsigpyip=sqrt(T(3,3));
betayip=p(4);dbetayip=dp(4);
alphyip=p(5);dalphyip=dp(5);

% Print data to screen and plot if requested
if printData
  disp(txt')
  figure(1)
  plot(S(id),1e6*sigxf,'b--')
  hold on
  plot_barsc(S(ido)',1e6*sigx',1e6*dsigx','b','o')
  hold off
  set(gca,'XLim',[S(id1),S(id2)])
  title('EXT Diagnostics Section')
  ylabel('Horizontal Beam Size (um)')
  xlabel('S (m)')
  plot_magnets_Lucretia(BEAMLINE(id),1,1);
  figure(2)
  plot(S(id),1e6*sigyf,'b--')
  hold on
  plot_barsc(S(ido)',1e6*sigy',1e6*dsigy','b','o')
  hold off
  set(gca,'XLim',[S(id1),S(id2)])
  title('EXT Diagnostics Section')
  ylabel('Vertical Beam Size (um)')
  xlabel('S (m)')
  plot_magnets_Lucretia(BEAMLINE(id),1,1);
else
  emitData=[energy ...
    emitx demitx emitxn demitxn embmxn dembmxn ...
    bmagx dbmagx bcosx dbcosx bsinx dbsinx ...
    betax dbetax bx0 alphx dalphx ax0 chi2x ...
    emity demity emityn demityn embmyn dembmyn ...
    bmagy dbmagy bcosy dbcosy bsiny dbsiny ...
    betay dbetay by0 alphy dalphy ay0 chi2y ...
    ido length(id) id' length(S) S' sigxf' sigx dsigx sigyf' sigy dsigy xf' yf'];
  for n=1:notr
    emitData=[emitData reshape(R{n},1,[])]; % R{n}=reshape(...,4,[])
  end
  emitData=[emitData ...
    exip0 bxip0 axip0 eyip0 byip0 ayip0 ...
    sigxip dsigxip sigpxip dsigpxip betaxip dbetaxip alphxip dalphxip ...
    sigyip dsigyip sigpyip dsigpyip betayip dbetayip alphyip dalphyip];

  if dointrinsic
    save(sprintf('userData/emit2dOTR_%s',datestr(now,30)), ...
      'rawotrdata','emitData','otruse','BEAMLINE','PS', ...
      'DX','DPX','DY','DPY','dDX','dDPX','dDY','dDPY', ...
      'ictdata','icterrdata');
  end
end
