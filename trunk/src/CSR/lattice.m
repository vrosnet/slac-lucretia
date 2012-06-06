function [Beam bo]=lattice
global BEAMLINE


% ==========
% Setup simple SBEN + DRIFT lattice
% ==========
BEAMLINE={};
Angle=2*asin(0.5/(2*1.5));
B = SBendStruc( 0.5, 2, Angle, [0 0], [0 0], [0.0064 0.0064], 0, 0, 'BEN1' );
D = DrifStruc( 0.5, 'DRIF1' );
BEAMLINE{1}=B; BEAMLINE{2}=D;
BEAMLINE{1}.P=0.15;
BEAMLINE{2}.P=0.15;
BEAMLINE{2}.S=B.L;
SetSPositions( 1, length(BEAMLINE), 0 );
BEAMLINE{1}.TrackFlag.LorentzDelay=1;
BEAMLINE{2}.TrackFlag.LorentzDelay=1;

% ============
% Make a beam structure
% ============
I=InitCondStruc;
I.Momentum=0.15;
I.x.NEmit=1e-99;
I.x.Twiss.beta=1;
I.x.Twiss.alpha=0;
I.y.NEmit=1e-99;
I.y.Twiss.beta=1;
I.y.Twiss.alpha=0;
I.sigz=50e-6;
I.SigPUncorrel=0;
I.PZCorrel=0;
I.Q=1e-9;
Beam=MakeBeam6DGauss(I,1e5,5,0);
Beam.Bunch.x(1:4,:)=0;

% =========
% Split BEAMLINE elements up and flag for CSR treatment
%   eleSplit(# segments to split SBEN into, # segments to split up
%   downstream CSR treatment area into, RMS bunch length to use for
%   estimation of length of downstream CSR region)
% =========
% eleSplitCSR(50,50,I.sigz);

% Set CSR flags Split, CSR_DriftSplit
for ibl=1:length(BEAMLINE)
  if strcmp(BEAMLINE{ibl}.Class,'SBEN')
    BEAMLINE{ibl}.TrackFlag.CSR=600;
    BEAMLINE{ibl}.TrackFlag.CSR_SmoothFactor=3;
    BEAMLINE{ibl}.TrackFlag.CSR_DriftSplit=50;
    BEAMLINE{ibl}.TrackFlag.Split=2;
    BEAMLINE{ibl}.TrackFlag.SynRad=2;
  end
end

% Beam tracking (NO CSR)
[stat bo]=TrackThru(1,length(BEAMLINE),Beam,1,1,0);
disp('Initial Beam:')
fprintf('sigma_z: %g sigma_E: %g mean E: %g\n',std(Beam.Bunch.x(5,:)),std(Beam.Bunch.x(6,:)),mean(Beam.Bunch.x(6,:)))
disp('CSR on:')
fprintf('sigma_z: %g sigma_E: %g mean E: %g\n',std(bo.Bunch.x(5,:)),std(bo.Bunch.x(6,:)),mean(bo.Bunch.x(6,:)))
