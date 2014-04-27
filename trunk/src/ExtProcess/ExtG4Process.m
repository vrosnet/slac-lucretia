classdef ExtG4Process < ExtProcess & ExtGeometry & handle
  %EXTG4PROCESS - class to handle interface of Lucretia beam using GEANT4 engine
  
  properties
    Verbose=0; % Verbosity level for printing status info to stdout in GEANT
    Ecut=0; % Energy cut for GEANT generated tracks / GeV
  end
  properties(Access=private)
    apercheck=true;
  end
  methods
    function obj = ExtG4Process(varargin)
      % Superclass initialization
      obj = obj@ExtProcess() ;
      obj = obj@ExtGeometry() ;
      % Get list of materials from GEANT4 database
      dbfile=which('G4MaterialsDatabase.txt');
      if isempty(dbfile)
        error('G4MaterialsDatabase.txt database file not on Matlab search path')
      end
      fid=fopen(dbfile,'r');
      while 1
        tline=fgetl(fid);
        if ~ischar(tline); break; end;
        t=regexp(tline,'(G4_\S+)','tokens','once');
        if ~isempty(t)
          obj.allowedMaterials{end+1}=t{1};
        end
      end
      fclose(fid);
      % Set other requested properties
      if ~nargin; return; end;
      if mod(nargin,2)
        error('Must supply property, value pairs as creation arguments')
      end
      for iarg=1:2:nargin
        if isprop(ExtG4Process,varargin{iarg})
          obj.(varargin{iarg})=varargin{iarg+1};
        else
          error('No settable property: %s',varargin{iarg})
        end
      end
      % Check any set properties
      obj.apercheck=false;
      checkExtGeometryProps(obj);
      obj.apercheck=true;
    end
  end
  methods
    function SetMaterial(obj,material)
      SetMaterial@ExtGeometry(obj,material);
      obj.Material=obj.Material;
    end
    function SetGeometryType(obj,type)
      global BEAMLINE
       if ~exist('type','var') || isempty(type)
         if ~isempty(obj.elemno) && isfield(BEAMLINE{obj.elemno},'Geometry')
           type=BEAMLINE{obj.elemno}.Geometry;
         else
           type='Ellipse';
         end
       end
       SetGeometryType@ExtGeometry(obj,type);
    end
    function SetAper(obj,val1,val2)
      global BEAMLINE
      % If this object attached to a BEAMLINE element then issue a warning
      % if the aperture does not match the BEAMLINE element or set to
      % BEAMLINE aper values if no values provided
      if ~exist('val1','var') || isempty(val1) && ~isempty(obj.elemno) && isfield(BEAMLINE{obj.elemno},'aper')
        val1=BEAMLINE{obj.elemno}.aper(1);
      end
      if ~exist('val2','var') || isempty(val2) && ~isempty(obj.elemno) && isfield(BEAMLINE{obj.elemno},'aper') && ...
        length(BEAMLINE{obj.elemno}.aper)>1
        val2=val1;
      elseif ~isempty(obj.elemno) && isfield(BEAMLINE{elemno},'aper') && length(BEAMLINE{obj.elemno}.aper)==1
        val2=val1;
      end
      SetAper@ExtGeometry(obj,val1,val2);
      if isempty(obj.elemno) || ~obj.apercheck; return; end;
      if isfield(BEAMLINE{obj.elemno},'aper')
        if BEAMLINE{obj.elemno}.aper(1)~=val1 || (length(BEAMLINE{obj.elemno}.aper)>1 && BEAMLINE{obj.elemno}.aper(2)~=val2)
          warning('ExtG4Process:AperMismatch','Provided geometry apertures don''t match associated Lucretia BEAMLINE element')
        end
      end
    end
  end
end

