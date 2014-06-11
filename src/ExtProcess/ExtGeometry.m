classdef ExtGeometry < handle
  %EXTGEOMETRY - Class to descibe geometry required for operation of
  %external processes (e.g. for GEANT4)
  
  properties(Constant, Hidden)
    allowedGeometryTypes={'Rectangle','Ellipse','Tapered'};
  end
  properties(Access=protected)
    allowedMaterials={'Vacuum'};
  end
  properties(SetAccess=protected)
    GeometryType='Ellipse'; % type of shape (e.g. 'Rectangle')
    AperX=1; % Inside half-aperture of shape (m) - Horizontal dimension
    AperY=1; % Inside half-aperture of shape (m) - Vertical dimension
    AperX2=0; % Secondary aperture definition (m) (for Tapered type)
    AperY2=0; % Secondary aperture definition (m) (for Tapered type)
    AperX3=0; % Secondary aperture definition (m) (for Tapered type)
    AperY3=0; % Secondary aperture definition (m) (for Tapered type)
    CollDX=1; % Collimator half-width (m) (for Tapered type)
    CollDY=1; % Collimator half-height (m) (for Tapered type)
    CollLen2=0; % Aux collimator length (for Tapered type)
    Thickness=1; % Geometry thickness (m) - defines half-aperture of World box
    Material='Vacuum'; % material type
    Material2='Vacuum'; % secondary material type
  end
  
  methods
    function obj=ExtGeometry(varargin)
      % Check inputs
      if mod(nargin,2)
        error('Must supply property, value pairs as creation arguments')
      end
      % parse inputs
      for iarg=1:2:nargin
        if isprop(ExtGeometry,varargin{iarg})
          obj.(varargin{iarg})=varargin{iarg+1};
        else
          error('No such property: %s',varargin{iarg})
        end
      end
      try
        obj.checkExtGeometryProps();
      catch ME
        error('Error constructing geometry:\n%',ME.message)
      end
    end
    function SetAper(obj,val1,val2)
      if ~exist('val1','var') || ~exist('val2','var')
        error('Must provide 2 values for AperX and AperY')
      end
      if any([val1 val2]<0)
        error('Apertures must be >=0')
      end
      obj.AperX=val1;
      obj.AperY=val2;
    end
    function SetAper2(obj,val1,val2)
      if ~exist('val1','var') || ~exist('val2','var')
        error('Must provide 2 values for AperX2 and AperY2')
      end
      if any([val1 val2]<0)
        error('Apertures must be >=0')
      end
      if any([val1 val2]>1)
        error('Apertures must be <=1')
      end
      obj.AperX2=val1;
      obj.AperY2=val2;
    end
    function SetAper3(obj,val1,val2)
      if ~exist('val1','var') || ~exist('val2','var')
        error('Must provide 2 values for AperX3 and AperY3')
      end
      if any([val1 val2]<0)
        error('Apertures must be >=0')
      end
      if any([val1 val2]>1)
        error('Apertures must be <=1')
      end
      obj.AperX3=val1;
      obj.AperY3=val2;
    end
    function SetGeometryType(obj,type)
      if ~exist('type','var') || ~ischar(type)
        error('Must provide Geometry Type String')
      end
      if ~any(strcmp(type,obj.allowedGeometryTypes))
        disp('Allowed Types:')
        disp(obj.allowedGeometryTypes)
        error('Must choose from list of allowed Geometry types ^')
      end
      obj.GeometryType=type;
    end
    function SetMaterial(obj,material)
      if ~any(strcmp(material,obj.allowedMaterials))
        error('Material not found in database');
      end
      obj.Material=material;
    end
    function SetMaterial2(obj,material)
      if ~any(strcmp(material,obj.allowedMaterials))
        error('Material not found in database');
      end
      obj.Material2=material;
    end
    function SetThickness(obj,val)
      if ~exist('val','var') || val<=0
        error('Must provide thickness parameter >0 (m)')
      end
      try
        obj.SetCollDX(obj.CollDX);
      catch
        obj.CollDX=val;
        if strcmp(obj.GeometryType,'Tapered')
          display('WARNING: Setting thickness to < CollDX, changing CollDX=Thickness')
        end
      end
      try
        obj.SetCollDY(obj.CollDY);
      catch
        obj.CollDY=val;
        if strcmp(obj.GeometryType,'Tapered')
          display('WARNING: Setting thickness to < CollDY, changing CollDY=Thickness')
        end
      end
      obj.Thickness=val;
    end
    function SetCollDX(obj,val)
      if ~exist('val','var') || val<=0 || val>obj.Thickness
        error('Must provide DX parameter >0 & <= obj.Thickness (m)')
      end
      obj.CollDX=val;
    end
    function SetCollDY(obj,val)
      if ~exist('val','var') || val<=0 || val>obj.Thickness
        error('Must provide DY parameter >0 & <= obj.Thickness (m)')
      end
      obj.CollDY=val;
    end
    function checkExtGeometryProps(obj)
      obj.SetGeometryType(obj.GeometryType);
      obj.SetAper(obj.AperX,obj.AperY);
      obj.SetAper2(obj.AperX2,obj.AperY2);
      obj.SetAper3(obj.AperX3,obj.AperY3);
      obj.SetThickness(obj.Thickness);
      obj.SetMaterial(obj.Material);
      obj.SetMaterial2(obj.Material2);
      obj.SetCollDX(obj.CollDX);
      obj.SetCollDY(obj.CollDY);
    end
  end
  
end

