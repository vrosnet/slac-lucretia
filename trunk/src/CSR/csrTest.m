global BEAMLINE
% =============
% Load in test lattice (SBEN + DRIFT)
% =============
Beam=lattice;

% ============
% Generate Track object
% ============
T=Track(Beam);
T.verbose=1;
T.csrStoreData=true;

% ===========
% Number of bins to use when histogramming for CSR
% ===========
T.csrNbins=300;
T.csrSmoothVal='robust';

% Beam tracking
T.trackThru;
B=T.beamOut;
disp('Split lattice + CSR:')
fprintf('sigma_z: %g sigma_E: %g <E>: %g\n',...
  std(B.Bunch.x(5,:)),std(B.Bunch.x(6,:)),mean(B.Bunch.x(6,:)))
me=[]; de=[]; L=[];
if T.csrStoreData
  for im=1:length(T.csrData)
    B=T.csrData(im).beam;
    de(im)=std(B.Bunch.x(6,:));
    me(im)=mean(Beam.Bunch.x(6,:))-mean(B.Bunch.x(6,:));
    L(im)=BEAMLINE{T.csrData(im).index}.S+BEAMLINE{T.csrData(im).index}.L;
  end
  figure
  plot(L,me.*1e3,'b-*',L,de.*1e3,'r-*')
end