function GenerateKnobs
clear global GIRDER BEAMLINE PS VERB MDATA INSTR FL
global GIRDER BEAMLINE PS VERB MDATA INSTR FL %#ok<NUSED>
load latticeFiles/ATF2lat
load latticeFiles/archive/ATF2-09May15_051906.mat

% Get track points
strack=findcells(BEAMLINE,'Name','IEX');
etrack=findcells(BEAMLINE,'Name','MW1IP');
% etrack=findcells(BEAMLINE,'Name','IP');

bs_target=5e-6;

if isempty(etrack)
  etrack=length(BEAMLINE);
  fudgeip=true;
else
  fudgeip=false;
end


if ~isfield(FL,'SimModel'); FL.SimModel=Model; end; %#ok<NODEF>

for isext=1:5
  sgi(isext)=BEAMLINE{FL.SimModel.magGroup.Sext.dB.ClusterList(isext).index(1)}.Girder;
%   PS(BEAMLINE{FL.SimModel.magGroup.Sext.dB.ClusterList(isext).index(1)}.Girder).Ampl=1;
end

if fudgeip
  L=0.4; % m
  T_ip=[1 L 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 L 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1];
else
  T_ip=[];
end

[stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0); %#ok<NODEF>
if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
fprintf('Target beamsize (Initial): %g\n',std(B1.Bunch.x(3,:)))
sigma=cov(B1.Bunch.x');
R=diag(ones(1,6));L=zeros(6,6);L(1,2)=1;L(3,4)=1;
xwaist=fminsearch(@(x) minWaist(x,R,L,sigma,1),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
ywaist=fminsearch(@(x) minWaist(x,R,L,sigma,3),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
fprintf('Target waist offset (Initial): %g (x) %g (y)\n',xwaist,ywaist)

% Generate Sext moves -> IP params matrix
% Switch off things I'm tired of hearing about
warning off MATLAB:lscov:RankDefDesignMat
warning off MATLAB:PSTrim_online
warning off MATLAB:MoverTrim_online
warning off MATLAB:FlHwUpdate_nonFSsim
Model.Measure_ipbcor=false; % turn on IP beam correlation measurement
Model.genKnobs=true;
nmoves=5; move_range=100e-6; 
Smoves=0:move_range/(nmoves-1):move_range; 
maxSEXT=max(abs(Smoves));
T=zeros(5,10);
sig3_x=zeros(5,nmoves);
sig3_y=zeros(5,nmoves);
Model.ipLaserRes=[0 0];
tic;
movefac=[1 1 1 1 1;
         0.5 0.5 0.5 0.5 0.5];
sf6=findcells(BEAMLINE,'Name','SF6FF');
% make new beams from fitted IEX parameters bx 13.61535 ax 1.36861 by
% 1.45331 ay -1.08288
Model.Initial_IEX.y.NEmit=28e-12*(1.28/0.511e-3);
Model.Initial_IEX.y.Twiss.alpha=-1.08288;
Model.Initial_IEX.y.Twiss.beta=1.45331;
Model.Initial_IEX.x.Twiss.NEmit=1.3e-9*(1.28/0.511e-3);
Model.Initial_IEX.x.Twiss.alpha=1.36831;
Model.Initial_IEX.x.Twiss.beta=13.61535;
Beam1_IEX = MakeBeam6DGauss( Model.Initial_IEX, 10000, 5, 1 );
Beam0_IEX = Beam1_IEX;
Beam0_IEX.Bunch.x=mean(Beam1_IEX.Bunch.x,2);
Beam0_IEX.Bunch.Q=mean(Beam1_IEX.Bunch.Q,2);
Beam0_IEX.Bunch.stop=mean(Beam1_IEX.Bunch.stop,2);

% Find good initial orbit
xmin=fminsearch(@(x) minorbit(x,Beam0_IEX,strack,etrack),[0 0 0 0 0],optimset('Display','none','MaxFunEvals',20000,'MaxIter',4000));
fprintf('Initial trajectory at IEX: %g / %g / %g / %g / %g(x/x''/y/y''/dE)\n',xmin)
Beam1_IEX.Bunch.x(1,:)=Beam1_IEX.Bunch.x(1,:)+xmin(1); 
Beam1_IEX.Bunch.x(2,:)=Beam1_IEX.Bunch.x(2,:)+xmin(2);
Beam1_IEX.Bunch.x(3,:)=Beam1_IEX.Bunch.x(3,:)+xmin(3);
Beam1_IEX.Bunch.x(4,:)=Beam1_IEX.Bunch.x(4,:)+xmin(4);
Beam1_IEX.Bunch.x(6,:)=Beam1_IEX.Bunch.x(6,:)+xmin(5);
[stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
fprintf('Target beamsize (Initial): %g\n',std(B1.Bunch.x(3,:)))
sigma=cov(B1.Bunch.x');
R=diag(ones(1,6));L=zeros(6,6);L(1,2)=1;L(3,4)=1;
xwaist=fminsearch(@(x) minWaist(x,R,L,sigma,1),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
ywaist=fminsearch(@(x) minWaist(x,R,L,sigma,3),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
fprintf('Target waist offset (Initial): %g (x) %g (y)\n',xwaist,ywaist)


% Put skew quad in kex2
kex2b=findcells(BEAMLINE,'Name','KEX2*'); kex=kex2b(end)-1;
qk4=findcells(BEAMLINE,'Name','QK4X'); qk4=qk4(end);
kex_s=BEAMLINE{kex}.S; kex_name=BEAMLINE{kex}.Name;
BEAMLINE{kex}=BEAMLINE{qk4};
BEAMLINE{kex}.Name=kex_name;
BEAMLINE{kex}.Block=[kex2b(1) kex2b(2)];
BEAMLINE{kex}.S=kex_s;
BEAMLINE{kex}=rmfield(BEAMLINE{kex},'Girder');
BEAMLINE{kex}=rmfield(BEAMLINE{kex},'PS');
AssignToPS( 1376, length(PS)+1 );
PS(end).Ampl=0; PS(end).SetPt=0;

% centre on beam
% [stat,B1]=TrackThru(Model.extStart,1375,Beam1_IEX,1,1,0);
% xpos=mean(B1.Bunch.x(1,:)); ypos=mean(B1.Bunch.x(3,:));
% BEAMLINE{1376}.Offset(1)=xpos;
% BEAMLINE{1376}.Offset(3)=ypos;

% Find waist
% qd0=findcells(BEAMLINE,'Name','QD0FF');
% qd0ps=BEAMLINE{qd0(1)}.PS;
% minqd0=fminbnd(@(x) minqd0(x,Beam1_IEX,strack,etrack,qd0ps,T_ip),0.9,1.1,optimset('Display','iter','MaxFunEvals',1000,'MaxIter',500));
% fprintf('Min QD0 PS: %g\n',minqd0)
% PS(qd0ps).Ampl=minqd0;

% Fit for incoming dispersion & coupling
% [stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
% if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
% fprintf('Target beamsize (pre disp/coupling optimisation): %g\n',std(B1.Bunch.x(3,:)))
% xmin=fminsearch(@(x) dispcoupfit(x,strack,etrack,T_ip,Model.Initial_IEX,bs_target,length(PS)),[0 0 0],optimset('Display','iter','MaxFunEvals',1000,'MaxIter',500,'TolFun',1e-2,'TolX',1e-5));
% Model.Initial_IEX.y.Twiss.eta=xmin(1);
% Model.Initial_IEX.y.Twiss.etap=xmin(2);
% Beam1_IEX = MakeBeam6DGauss( Model.Initial_IEX, 10000, 5, 1 );
% Beam0_IEX = Beam1_IEX;
% Beam0_IEX.Bunch.x=mean(Beam1_IEX.Bunch.x,2);
% Beam0_IEX.Bunch.Q=mean(Beam1_IEX.Bunch.Q,2);
% Beam0_IEX.Bunch.stop=mean(Beam1_IEX.Bunch.stop,2);
% PS(end).Ampl=xmin(3);
% [stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
% if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
% fprintf('Target beamsize (post disp/coupling optimisation): %g\n',std(B1.Bunch.x(3,:)))
% sigma=cov(B1.Bunch.x');
% R=diag(ones(1,6));L=zeros(6,6);L(1,2)=1;L(3,4)=1;
% xwaist=fminsearch(@(x) minWaist(x,R,L,sigma,1),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
% ywaist=fminsearch(@(x) minWaist(x,R,L,sigma,3),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
% fprintf('Target waist offset: %g (x) %g (y)\n',xwaist,ywaist)
% qk4=findcells(BEAMLINE,'Name','QK1X');
% qk4ps=BEAMLINE{qk4(1)}.PS;
% minqk4=fminbnd(@(x) minsq(x,Beam1_IEX,strack,etrack,qk4ps,T_ip,bs_target),-1,1,optimset('Display','iter','MaxFunEvals',1000,'MaxIter',500,'TolFun',1e-5,'TolX',1e-5));
% PS(qk4ps).Ampl=minqk4;
% [stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
% if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
% fprintf('Target beamsize (post qk1 min): %g\n',std(B1.Bunch.x(3,:)))
% qk4=findcells(BEAMLINE,'Name','QK4X');
% qk4ps=BEAMLINE{qk4(1)}.PS;
% minqk4=fminbnd(@(x) minsq(x,Beam1_IEX,strack,etrack,qk4ps,T_ip,bs_target),-1,1,optimset('Display','iter','MaxFunEvals',1000,'MaxIter',500,'TolFun',1e-5,'TolX',1e-5));
% PS(qk4ps).Ampl=minqk4;
% [stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
% if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
% fprintf('Target beamsize (post qk4 min): %g\n',std(B1.Bunch.x(3,:)))
% qk4=findcells(BEAMLINE,'Name','QK2X');
% qk4ps=BEAMLINE{qk4(1)}.PS;
% minqk4=fminbnd(@(x) minsq(x,Beam1_IEX,strack,etrack,qk4ps,T_ip,bs_target),-1,1,optimset('Display','iter','MaxFunEvals',1000,'MaxIter',500,'TolFun',1e-5,'TolX',1e-5));
% PS(qk4ps).Ampl=minqk4;
% [stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
% if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
% fprintf('Target beamsize (post qk2 min): %g\n',std(B1.Bunch.x(3,:)))
% qk4=findcells(BEAMLINE,'Name','QK3X');
% qk4ps=BEAMLINE{qk4(1)}.PS;
% minqk4=fminbnd(@(x) minsq(x,Beam1_IEX,strack,etrack,qk4ps,T_ip,bs_target),-1,1,optimset('Display','iter','MaxFunEvals',1000,'MaxIter',500,'TolFun',1e-5,'TolX',1e-5));
% PS(qk4ps).Ampl=minqk4;
% [stat,B1]=TrackThru(strack,etrack,Beam1_IEX,1,1,0);
% if fudgeip; B1.Bunch.x=T_ip*B1.Bunch.x; end;
% fprintf('Target beamsize (post qk3 min): %g\n',std(B1.Bunch.x(3,:)))

% make knobs
[stat,B1]=TrackThru(strack,sf6(1)-1,Beam1_IEX,1,1,0);
for iSext=1:5
  GIRDER{sgi(iSext)}.MoverPos=[0 0 0];
  fprintf('Generating knobs for Sext %d (x)\n',iSext)
  % move in x
  for imove=1:length(Smoves)
    GIRDER{sgi(iSext)}.MoverPos=[Smoves(imove)*movefac(1,iSext) 0 0];
    % Get IP params
    [stat,B]=TrackThru(sf6(1),etrack,B1,1,1,0);
    if fudgeip; sigma=T_ip*cov(B.Bunch.x')*T_ip'; else sigma=cov(B.Bunch.x'); end;
    sig3_x(iSext,imove) = sigma(2,3);
    R=diag(ones(1,6));L=zeros(6,6);L(1,2)=1;L(3,4)=1;
    xWaist_x(iSext,imove)=fminsearch(@(x) minWaist(x,R,L,sigma,1),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
    xWaist_y(iSext,imove)=fminsearch(@(x) minWaist(x,R,L,sigma,3),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
    [xDisp_x(iSext,imove) xDisp_y(iSext,imove)] = getdisp(B1,sf6(1),etrack,T_ip);
  end
  % Move in y
  fprintf('Generating knobs for Sext %d (y)\n',iSext)
  for imove=1:length(Smoves)
    GIRDER{sgi(iSext)}.MoverPos=[0 Smoves(imove)*movefac(1,iSext) 0];
    % Get IP params
    [stat,B]=TrackThru(sf6(1),INSTR{end}.Index,B1,1,1,0);
    if fudgeip; sigma=T_ip*cov(B.Bunch.x')*T_ip'; else sigma=cov(B.Bunch.x'); end;
    sig3_y(iSext,imove) = sigma(2,3);
    R=diag(ones(1,6));L=zeros(6,6);L(1,2)=1;L(3,4)=1;
    yWaist_x(iSext,imove)=fminsearch(@(x) minWaist(x,R,L,sigma,1),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
    yWaist_y(iSext,imove)=fminsearch(@(x) minWaist(x,R,L,sigma,3),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
    [yDisp_x(iSext,imove) yDisp_y(iSext,imove)] = getdisp(B1,sf6(1),etrack,T_ip);
  end
end
GIRDER{sgi(iSext)}.MoverPos=[0 0 0];
%   Linear fits
for iSext=1:5
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(1,1),xWaist_x(iSext,:)./0.004,1,1); T(1,(iSext-1)*2+1)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(2,1),yWaist_x(iSext,:)./0.004,1,1); T(1,(iSext-1)*2+2)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(1,2),xWaist_y(iSext,:)./1e-4,1,1); T(2,(iSext-1)*2+1)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(2,2),yWaist_y(iSext,:)./1e-4,1,1); T(2,(iSext-1)*2+2)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(1,3),xDisp_x(iSext,:)./1e-5,1,1); T(3,(iSext-1)*2+1)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(2,3),yDisp_x(iSext,:)./1e-5,1,1); T(3,(iSext-1)*2+2)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(1,4),xDisp_y(iSext,:)./1e-5,1,1); T(4,(iSext-1)*2+1)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(2,4),yDisp_y(iSext,:)./1e-5,1,1); T(4,(iSext-1)*2+2)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(1,5),sig3_x(iSext,:)./5e-11,1,1); T(5,(iSext-1)*2+1)=q(2);
  q=noplot_polyfit((Smoves./maxSEXT).*movefac(2,5),sig3_y(iSext,:)./5e-11,1,1); T(5,(iSext-1)*2+2)=q(2);
end % for iSext

% Make knobs
Knobs.WaistX.smoves = lscov(T,[1 0 0 0 0]').*maxSEXT;
Knobs.WaistX.scale = 0.004;
Knobs.WaistY.smoves = lscov(T,[0 1 0 0 0]').*maxSEXT;
Knobs.WaistY.scale = 1e-4;
Knobs.DispX.smoves = lscov(T,[0 0 1 0 0]').*maxSEXT;
Knobs.DispX.scale = 1e-5;
Knobs.DispY.smoves = lscov(T,[0 0 0 1 0]').*maxSEXT;
Knobs.DispY.scale = 1e-5;
Knobs.XpY.smoves = lscov(T,[0 0 0 0 1]').*maxSEXT;
Knobs.XpY.scale = 5e-11;
save sextknobs Knobs
%%
% test
ik1{1}=reshape(Knobs.DispY.smoves./max(abs(Knobs.DispY.smoves)),2,5);
ik1{2}=reshape(Knobs.XpY.smoves./max(abs(Knobs.XpY.smoves)),2,5);
ik1{3}=reshape(Knobs.WaistY.smoves./max(abs(Knobs.WaistY.smoves)),2,5);
k=-1500:500:1500;
knames={'Dispersion' 'Coupling' 'Waist'};
for iknob=1:3
  fprintf('Testing Knob: %s\n',knames{iknob})
  for ik=1:length(k)
    for isext=1:5
      GIRDER{sgi(isext)}.MoverPos=[ik1{iknob}(1,isext) ik1{iknob}(2,isext) 0].*k(ik).*1e-6;
    end
    [stat,B]=TrackThru(sf6(1),etrack,B1,1,1,0);
    if fudgeip; B.Bunch.x=T_ip*B.Bunch.x; end;
    sigma=cov(B.Bunch.x');
    sig3(iknob,ik) = sigma(2,3);
    ysize(iknob,ik)=std(B.Bunch.x(3,:)); 
    R=diag(ones(1,6));L=zeros(6,6);L(1,2)=1;L(3,4)=1;
    xwaist(iknob,ik)=fminsearch(@(x) minWaist(x,R,L,sigma,1),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
    ywaist(iknob,ik)=fminsearch(@(x) minWaist(x,R,L,sigma,3),0,optimset('Tolx',1e-6,'TolFun',0.1e-6^2));
    [dx dy] = getdisp(B1,sf6(1),etrack,T_ip); dispx(iknob,ik)=dx; dispy(iknob,ik)=dy;
  end
end
figure; plot(k,ysize')
figure
subplot(3,2,1), plot(dispx'), title('dispx')
subplot(3,2,2), plot(dispy'), title('dispy')
subplot(3,2,4), plot(sig3'), title('sig23')
subplot(3,2,5), plot(xwaist'), title('xwaist')
subplot(3,2,6), plot(ywaist'), title('ywaist')
save

function chi2 = minWaist(x,R,L,sig,dir)

newsig=(R+L.*x(1))*sig*(R+L.*x(1))';
chi2=newsig(dir,dir)^2;

function chi2 = sinFit(x,data,error) %#ok<DEFNU>

chi2=sum( ( data - ( x(1) * sin((1:length(data))/x(2)+2*pi*x(3))+mean(data) ) ).^2 ./ error.^2);



function chi2=minqd0(x,beam,istart,iend,iqd0,T) %#ok<DEFNU>
global PS

PS(iqd0).Ampl=x;
[stat,B]=TrackThru(istart,iend,beam,1,1,0);
if ~isempty(T); B.Bunch.x=T*B.Bunch.x; end;
chi2=std(B.Bunch.x(3,:))^2/1e-6^2;

function chi2=minsq(x,beam,istart,iend,ind,T,bs_target)
global PS
persistent xmin

if isempty(xmin)
  xmin=[0 0 0 0 0];
end

PS(ind).Ampl=x;
[stat,B]=TrackThru(istart,iend,beam,1,1,0);
if ~isempty(T); B.Bunch.x=T*B.Bunch.x; end;
chi2=abs((std(B.Bunch.x(3,:))-bs_target)^2/1e-6^2);

function chi2 = mintwiss(x,Initial,istart,iend,T)
Initial.x.Twiss.alpha=x(1);
Initial.x.Twiss.beta=x(2);
Initial.y.Twiss.alpha=x(3);
Initial.y.Twiss.beta=x(4);
beam = MakeBeam6DGauss( Initial, 10001, 3, 1 );
[stat,B]=TrackThru(istart,iend,beam,1,1,0);
B.Bunch.x=T*B.Bunch.x;
chi2=std(B.Bunch.x(3,:))^2/1e-6^2;

function chi2 = minorbit(x,beam,istart,iend)

beam.Bunch.x(1,:)=beam.Bunch.x(1,:)+x(1);
beam.Bunch.x(2,:)=beam.Bunch.x(2,:)+x(2);
beam.Bunch.x(3,:)=beam.Bunch.x(3,:)+x(3);
beam.Bunch.x(4,:)=beam.Bunch.x(4,:)+x(4);
beam.Bunch.x(6,:)=beam.Bunch.x(6,:)+x(5);

[stat,B,instdata]=TrackThru(istart,iend,beam,1,1,0);
xdat=[instdata{1}.x]; ydat=[instdata{1}.y];
if length(xdat)<34
  xdat(end+1:34)=5e-3;
  ydat(end+1:34)=5e-3;
end
chi2=sum(xdat.^2+ydat.^2)./sum(ones(1,34).*1e-6.^2);

function chi2 = dispcoupfit(x,istart,iend,T,Initial,bs_target,kexind)
global PS
PS(kexind).Ampl=x(3);
Initial.y.Twiss.eta=x(1);
Initial.y.Twiss.etap=x(2);
beam = MakeBeam6DGauss( Initial, 10001, 3, 1 );
[stat,B]=TrackThru(istart,iend,beam,1,1,0);
if ~isempty(T); B.Bunch.x=T*B.Bunch.x; end;
chi2=abs((std(B.Bunch.x(3,:))-bs_target)^2/1e-6^2);

function [xDisp yDisp] = getdisp(B1,ind1,ind2,T_ip)

% Get x and y IP dispersion
E_nominal=mean(B1.Bunch.x(6,:));
espread=0.5e-2;
ndisp=10;
disp=-espread/2:espread/(ndisp-1):espread/2;
beam1=B1;
ip_x=zeros(ndisp,1); ip_y=zeros(ndisp,1);
for idisp=1:ndisp
  beam1.Bunch.x(6,:)=B1.Bunch.x(6,:)+(E_nominal*disp(idisp));
  [stat,B]=TrackThru(ind1,ind2,beam1,1,1,0);
  if ~isempty(T_ip); B.Bunch.x=T_ip*B.Bunch.x; end;
  ip_x(idisp)=mean(B.Bunch.x(1,:)); ip_y(idisp)=mean(B.Bunch.x(3,:));
end
q=noplot_polyfit(disp,ip_x,1,1);
xDisp = q(2);
q=noplot_polyfit(disp,ip_y,1,1);
yDisp = q(2);