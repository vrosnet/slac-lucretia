classdef ExtSecondaryParticles
  %EXTSECONDARYPARTICLES - container for secondary particles generated by
  %external processes
  
  properties
    ParticleType % Cell array of strings describing particle types
    Pos % position/Energy vectors of particles (in Beam.Bunch.x format)
    PrimaryID % Link to primary macro particle in tracked Lucretia Beam
    ProcType % Link to physics process which generated secondary (list of processes in supportedProcessPhysics property)
    TrackStatus % Status of secondary tracks when written out (proccess dependent, for GEANT: listed status options in G4TrackStatus property)
    NumStored % Number of secondary particles for this bunch that were stored
  end
  
  methods
    function obj = ExtSecondaryParticles(maxpart)
      if ~exist('maxpart','var')
        maxpart=0;
      end
      obj.ParticleType=cell(1,maxpart);
      obj.Pos=zeros(6,maxpart);
      obj.PrimaryID=uint32(zeros(1,maxpart));
    end
  end
  
end

