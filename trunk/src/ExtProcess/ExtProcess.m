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
    SecondaryParticleTypes = 'All' ; % Which particle types to store
    PrimarySampleOrder
    PrimaryParticlesData
    SecondaryParticlesData
  end
  properties(Dependent)
    MaxSecondaryParticlesPerPrimary
    MaxSecondaryParticles % set >0 to store up to N secondary particles produced
    MaxPrimaryParticles
  end
  properties(Access=protected)
    fPrimarySampleOrder ;
    fMaxSecondaryParticlesPerPrimary = uint32(10) ;
    fMaxSecondaryParticles= uint32(0) ;
    fMaxPrimaryParticles = uint32(1e4) ;
    fNumSecondariesStored = uint32(0) ; % Number of secondary particles actually stored
  end
  
  methods
    function obj = ExtProcess()
    end
    function val=get.MaxPrimaryParticles(obj)
      val=obj.fMaxPrimaryParticles;
    end
    function set.MaxPrimaryParticles(obj,val)
      obj.fMaxPrimaryParticles=uint32(val);
    end
    function val=get.MaxSecondaryParticlesPerPrimary(obj)
      val=obj.fMaxSecondaryParticlesPerPrimary;
    end
    function val=get.MaxSecondaryParticles(obj)
      val=obj.fMaxSecondaryParticles;
    end
    function set.MaxSecondaryParticlesPerPrimary(obj,val)
      obj.fMaxSecondaryParticlesPerPrimary=uint32(val);
    end
    function set.MaxSecondaryParticles(obj,val)
      obj.fMaxSecondaryParticles=uint32(val);
    end
    function FinalizeTrackingData(obj)
      % To be called by Track after tracking
      global BEAMLINE
      if isfield(BEAMLINE{obj.elemno},'ExtProcess_primariesData')
        BEAMLINE{obj.elemno}.ExtProcess(obj.processID).PrimaryParticlesData=BEAMLINE{obj.elemno}.ExtProcess_primariesData;
        BEAMLINE{obj.elemno}=rmfield(BEAMLINE{obj.elemno},'ExtProcess_primariesData') ;
      end
      if isfield(BEAMLINE{obj.elemno},'ExtProcess_secondariesData')
        BEAMLINE{obj.elemno}.ExtProcess(obj.processID).SecondaryParticlesData=BEAMLINE{obj.elemno}.ExtProcess_secondariesData;
        BEAMLINE{obj.elemno}=rmfield(BEAMLINE{obj.elemno},'ExtProcess_secondariesData') ;
      end
    end
    function InitializeTrackingData(obj,primaryBeam,b1,b2)
      % Set order in which primary rays are serviced, order by highest
      % charge weight first, with equal charge weights randomized
      global BEAMLINE
      if ~isfield(primaryBeam,'Bunch') || ~isfield(primaryBeam.Bunch,'Q') || length(primaryBeam.Bunch.Q)<1
        error('Badly formatted Lucretia Beam ''primaryBeam''')
      end
      for ibunch=b1:b2
        randind=randperm(length(primaryBeam.Bunch(ibunch).Q));
        [~, sortind]=sort(primaryBeam.Bunch(ibunch).Q(randind));
        obj.PrimarySampleOrder{ibunch}=uint32(randind(sortind));
      end
      % Set up data structures for returning info about primary and
      % secondary particles not codified in Lucretia Beam
      for ibunch=b1:b2
        BEAMLINE{obj.elemno}.ExtProcess_primariesData(ibunch)=ExtPrimaryParticles();
        BEAMLINE{obj.elemno}.ExtProcess_secondariesData(ibunch)=ExtSecondaryParticles();
      end
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

