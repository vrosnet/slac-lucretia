classdef ExtPrimaryParticles
  %EXTPRIMARYPARTICLES - Data associated with primary particles handed over
  %to EXT processes
  
  properties
    regeneratedID % List of particles stopped by Lucretia tracking that EXT process re-started
  end
  
  methods
    function obj = ExtPrimaryParticles(maxpart)
      if ~exist('maxpart','var')
        maxpart=0;
      end
      obj.regeneratedID=uint32(zeros(1,maxpart));
    end
  end
  
end

