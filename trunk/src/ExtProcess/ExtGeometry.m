classdef ExtGeometry < handle
  %EXTGEOMETRY - Class to descibe geometry required for operation of
  %external processes (e.g. for GEANT4)
  
  properties(Constant)
    allowedGeometryTypes={'Rectangle','Ellipse'};
  end
  properties(SetAccess=protected)
    allowedMaterials={'Vacuum'};
    GeometryType='Ellipse'; % type of shape (e.g. 'Rectangule')
    AperX=1; % Inside half-aperture of shape (m) - Horizontal dimension
    AperY=1; % Inside half-aperture of shape (m) - Vertical dimension
    Thickness=1; % Geometry thickness (m)
    Material='Vacuum'; % material type
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
    function SetThickness(obj,val)
      if ~exist('val','var') || val<=0
        error('Must provide thickness parameter >0 (m)')
      end
      obj.Thickness=val;
    end
    function checkExtGeometryProps(obj)
      obj.SetGeometryType(obj.GeometryType);
      obj.SetAper(obj.AperX,obj.AperY);
      obj.SetThickness(obj.Thickness);
      obj.SetMaterial(obj.Material);
    end
  end
  
end

