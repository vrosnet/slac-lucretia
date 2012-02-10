% =============
% Load in test lattice (SBEN + DRIFT)
% =============
Beam=lattice;

% ============
% Generate Track object
% ============
T=Track(Beam);
T.verbose=0;
T.csrStoreData=true;

% ===========
% Number of bins to use when histogramming for CSR
% ===========
T.csrNbins=3000;

% Beam tracking
T.trackThru;
B=T.beamOut;
disp('Split lattice + CSR:')
fprintf('sigma_z: %g sigma_E: %g <E>: %g\n',...
  std(B.Bunch.x(5,:)),std(B.Bunch.x(6,:)),mean(B.Bunch.x(6,:)))


% function doplot(ibl,bsplit,bo,dE,me,me1,s,se,z)
% figure(1)
% if ibl<=bsplit
%   plot(-z./std(bo.Bunch.x(5,:)),dE.*1e3) % dE / MeV (in Bend)
% else
%   plot(-z./std(bo.Bunch.x(5,:)),dE.*1e3,'r') % dE / MeV (beyond Bend)
% end
% xlabel('s/\sigma_s'); ylabel('qW(s).s / MeV')
% hold on
% figure(2)
% plot(bo.Bunch.x(5,:),bo.Bunch.x(6,:).*1000,'.')
% xlabel('Z / m');ylabel('E / MeV')
% figure(3)
% plot(s,abs(me-me1).*1e3,'b',s,se.*1e3,'r')
% xlabel('(s-s_{bend}) / m'); ylabel('<\DeltaE> and \sigma_E / MeV')
% legend('<\DeltaE>','\sigma_E','Location','NorthWest')
% drawnow('expose')
% pause(0.5)
