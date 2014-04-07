classdef ExtProcess < handle
  % EXTPROCESS - Class for handling interfaces to external code linked with Lucretia
  % Currently supported is GEANT4
  
  properties(Constant, Hidden)
    supportedProcessTypes = {'GEANT4'} ;
    supportedProcessPhysics = {'All'} ;
    supportedParticleTypes = {'All'} ;
  end
  properties(Abstract)
    Verbose ;
  end
  properties(Access=protected)
    elemno ;
    processID ;
    type ;
    processPhysicsEnabled = 'All' ;
  end
  properties(SetAccess=protected)
    SecondaryBeam ; % Lucretia Beam structure for store secondaries
    SecondaryParticleTypes = 'All' ; % Which particle types to store
    PrimaryOrder % Order in which primary particles are serviced
  end
  properties
    MaxSecondaryParticles = 0 ; % set >0 to store up to N secondary particles produced
    MaxPrimaryParticles = 1e4 ;
    MaxSecondaryParticlesPerPrimary = 10 ;
    NumSecondariesStored = 0 ; % Number of secondary particles actually stored
  end
  
  methods
    function obj = ExtProcess()
    end
    function InitializeSecondaries(obj, primaryBeam)
      global BEAMLINE
      if obj.MaxSecondaryParticles<1; return; end;
      if ~isfield(primaryBeam,'Bunch') || ~isfield(primaryBeam.Bunch,'Q') || length(primaryBeam.Bunch.Q)<1
        error('Badly formatted Lucretia Beam ''primaryBeam''')
      end
      obj.SecondaryBeam = CreateBlankBeam(length(primaryBeam.Bunch),min([obj.MaxPrimaryParticles,obj.MaxSecondaryParticles]), ...
        BEAMLINE{obj.elemno}.P,primaryBeam.BunchInterval) ;
      obj.SecondaryBeam.Bunch.type = cell(1,length(obj.SecondaryBeam.Bunch.Q)) ;
    end
    function SetPrimaryOrdering(obj,primaryBeam)
      % Set order in which primary rays are serviced, order by highest
      % charge weight first, with equal charge weights randomized
      if ~isfield(primaryBeam,'Bunch') || ~isfield(primaryBeam.Bunch,'Q') || length(primaryBeam.Bunch.Q)<1
        error('Badly formatted Lucretia Beam ''primaryBeam''')
      end
      randind=randperm(length(primaryBeam.Bunch.Q));
      [~, sortind]=sort(primaryBeam.Bunch.Q(randind));
      obj.PrimaryOrder=randind(sortind);
    end
  end
  methods(Static)
    function new(processType,elemno,varargin)
      global BEAMLINE
      % First parse inputs and generate new ExtProcess structure in
      % BEAMLINE
      if elemno<1 || elemno>length(BEAMLINE)
        error('No existing BEAMLINE element matching request')
      end
      if ~isfield(BEAMLINE{elemno},'ExtProcess') || isempty(BEAMLINE{elemno}.ExtProcess)
        procno=1;
      else
        procno=length(BEAMLINE{elemno}.ExtProcess);
      end
      switch processType
        case 'GEANT4'
          % Check there is some length to this element
          if ~isfield(BEAMLINE{elemno},'L')
            error('No length field to this element, GEANT4 tracking not possible!')
          end
          if BEAMLINE{elemno}.L==0
            warning('ExtProcess:noLength','Zero length for this element, GEANT4 will not attempt tracking here');
          end
          BEAMLINE{elemno}.ExtProcess(procno)=ExtG4Process(varargin{:});
          BEAMLINE{elemno}.ExtProcess(procno).elemno=elemno;
          BEAMLINE{elemno}.ExtProcess(procno).type='GEANT4';
          % If not explicitly setting apertures, set them to default values
          % (equal to BEAMLINE aper field if there is one)
          if ~any(strcmp('AperX',varargin)) && ~any(strcmp('AperY',varargin))
            BEAMLINE{elemno}.ExtProcess(procno).SetAper() ;
          end
          % If not explicitly set and the BEAMLINE element has a Geometry
          % flag then set the GeometryType as this
          if ~any(strcmp('GeometryType',varargin))
            BEAMLINE{elemno}.ExtProcess(procno).SetGeometryType();
          end
        otherwise
          error('Unsupported EXT Process Type')
      end
      BEAMLINE{elemno}.ExtProcess(procno).processID=procno;
    end
  end
  
end

