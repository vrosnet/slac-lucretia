classdef physConsts < handle
  % Definition of physics constants
  
  properties(Constant)
    emass=0.51099906e-3; % electron rest mass (GeV)
    eQ=1.602176462e-19;  % electron charge (C)
    clight=2.99792458e8; % speed of light (m/sec)
    Cb=1e9/2.99792458e8; % rigidity conversion (T-m/GeV)
  end
end

