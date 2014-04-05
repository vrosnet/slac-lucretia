classdef ExtProcess < handle
  % EXTPROCESS - Class for handling interfaces to external code linked with Lucretia
  % Currently supported is GEANT4
  
  properties(Constant)
    supportedProcessTypes = {'GEANT4'} ;
    supportedProcessPhysics = {'All'} ;
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
  properties
    MaxPrimaryParticles = 1e4 ;
    MaxSecondaryParticlesPerPrimary = 10 ;
  end
  
  methods
    function obj = ExtProcess()
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

