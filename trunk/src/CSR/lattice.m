function Beam=lattice
global BEAMLINE


% ==========
% Setup simple SBEN + DRIFT lattice
% ==========
BEAMLINE={};
B = SBendStruc( 0.9125*2, 1.2586*2, 0.0164*2, [0.0164 0.0164], [0 0], [0.0064 0.0064], 0.5, 0, 'BEN1' );
D = DrifStruc( 6, 'DRIF1' );
BEAMLINE{1}=B; BEAMLINE{2}=D;
BEAMLINE{1}.P=23;
BEAMLINE{2}.P=23;
BEAMLINE{2}.S=B.L;
SetSPositions( 1, length(BEAMLINE), 0 );

% ============
% Make a beam structure
% ============
I=InitCondStruc;
I.Momentum=23;
I.x.NEmit=1e-30;
I.x.Twiss.beta=1;
I.x.Twiss.alpha=0;
I.y.NEmit=1e-30;
I.y.Twiss.beta=1;
I.y.Twiss.alpha=0;
I.sigz=51e-6;
I.SigPUncorrel=0;
I.PZCorrel=0;
I.Q=1.6022e-19*2e10;
Beam=MakeBeam6DGauss(I,1e5,5,0);

% =========
% Split BEAMLINE elements up and flag for CSR treatment
%   eleSplit(# segments to split SBEN into, # segments to split up
%   downstream CSR treatment area into, RMS bunch length to use for
%   estimation of length of downstream CSR region)
% =========
eleSplitCSR(50,50,I.sigz);

% Beam tracking (NO CSR)
[stat bo]=TrackThru(1,length(BEAMLINE),Beam,1,1,0);
disp('Initial Beam:')
fprintf('sigma_z: %g sigma_E: %g mean E: %g\n',std(Beam.Bunch.x(5,:)),std(Beam.Bunch.x(6,:)),mean(Beam.Bunch.x(6,:)))
disp('CSR off:')
fprintf('sigma_z: %g sigma_E: %g mean E: %g\n',std(bo.Bunch.x(5,:)),std(bo.Bunch.x(6,:)),mean(bo.Bunch.x(6,:)))