classdef ExtSecondaryParticles < handle
  %EXTSECONDARYPARTICLES - container for secondary particles generated by
  %external processes
  
  properties(SetAccess=protected)
    secondaryParticleType % Cell array of strings describing particle types
    secondaryPosVector % position/Energy vectors of particles (in Beam.Bunch.x format)
  end
  
  methods
  end
  
end

