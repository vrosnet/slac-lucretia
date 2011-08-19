classdef FlGui < handle
  %FLGUI Template for Floodland GUIs
  
  properties
    guiFont
    guiNumScreens
    guiScreenSize
  end
  properties
    gui
  end
  properties(Access=private)
    guiTitle
  end
  
  methods
    %% constructor
    function obj=FlGui
      mp=get(0,'MonitorPositions');
      sz=size(mp);
      obj.guiNumScreens=sz(1);
      obj.guiScreenSize=mp;
      obj.guiFont=get(0,'FixedWidthFontName');
    end
    %% test for gui existence
    function resp=guiExists(obj,name)
      resp = isprop(obj,'gui') && isfield(obj.gui,name) && ~isempty(obj.gui.(name)) && ishandle(obj.gui.(name)) && ...
        strcmp(get(obj.gui.(name),'Name'),obj.guiTitle.(name));
    end
    %% figure
    function guiCreateFigure(obj,name,title,size)
      if ~exist('size','var') || (length(size)~=2 && length(size)~=4)
        error('Size variable should be supplied as vector length 2 or 4');
      end
      if ~exist('name','var') || ~ischar(name)
        error('Must supply name for figure handle')
      end
      if ~exist('title','var') || ~ischar(name)
        error('Must supply title for figure')
      end
      % If GUI already exists, re-display it
      if obj.guiExists(name)
        oldpos=get(obj.gui.(name),'Position');
        close(obj.gui.(name))
        obj.gui.(name)=obj.(name);
        set(obj.gui.(name),'Position',oldpos)
        return
      end
      % if size is 1*2, default to making figure window in centre of
      % current screen
      % - find which screen we are in
      thisScreen=0;
      for iscreen=1:obj.guiNumScreens
        curpos=get(0,'PointerLocation');
        if inpolygon(curpos(1),curpos(2),[obj.guiScreenSize(iscreen,1) obj.guiScreenSize(iscreen,1)+obj.guiScreenSize(iscreen,3)],...
            [obj.guiScreenSize(iscreen,2) obj.guiScreenSize(iscreen,2)+obj.guiScreenSize(iscreen,4)])
          thisScreen=iscreen;
          break
        end
      end
      if ~thisScreen; error('Pointer outside screen area!'); end;
      figSize=[min([size(1) obj.guiScreenSize(thisScreen,3)]) min([(obj.guiScreenSize(thisScreen,4)) size(2)])];
      if length(size)==2 % put it in the middle of current screen
        xpos=obj.guiScreenSize(iscreen,1)+obj.guiScreenSize(thisScreen,3)/2;
        ypos=obj.guiScreenSize(iscreen,2)+obj.guiScreenSize(thisScreen,4)/2;
        pos=[xpos-figSize(1)/2 ypos+figSize(2)/2 figSize];
      else % otherwise put it exactly where requested
        pos=size;
      end
      obj.gui.(name)=figure('position',pos,'Name',title,'NumberTitle','off','MenuBar','none');
      obj.guiTitle.(name)=title;
    end
    %% panel
    function guiCreatePanel(obj,name,title,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name) = uipanel('Parent',parent,'FontName',obj.guiFont,'Title',title,'Units','normalized','Position',pos);
    end
    %% button group
    function guiCreateButtonGroup(obj,name,title,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name) = uibuttongroup('Parent',parent,'FontName',obj.guiFont,'Title',title,'Units','normalized','Position',pos);
    end
    %% axes
    function guiCreateAxes(obj,name,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=axes('Parent',parent,'Units','normalized','FontName',obj.guiFont,'Position',pos);
    end
    %% pushbutton
    function guiCreatePushbutton(obj,name,title,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','pushbutton','String',title,...
          'FontWeight','bold','FontName',obj.guiFont,'Units','normalized','Position',pos);
    end
    %% edit
    function guiCreateEdit(obj,name,txt,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','edit','FontName',obj.guiFont,'Units','normalized','String',txt,...
        'Position',pos,'BackgroundColor','white') ;
    end
    %% text
    function guiCreateText(obj,name,txt,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
       obj.gui.(name)=uicontrol(parent,'Style','text','FontName',obj.guiFont,'Units','normalized','String',txt,...
        'Position',pos,'FontWeight','bold') ;
    end
    %% popupmenu
    function guiCreatePopupmenu(obj,name,menutxt,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','popupmenu','String',menutxt,'FontName',obj.guiFont,'Units','normalized',...
        'Position',pos,'BackgroundColor','white','FontWeight','bold');
    end
    %% radiobutton
    function guiCreateRadiobutton(obj,name,txt,val,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','radiobutton','String',txt,'FontName',obj.guiFont,...
        'Units','normalized','Position',pos,'FontWeight','bold','Value',val);
    end
    %% listbox
    function guiCreateListbox(obj,name,txt,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','listbox','String',txt,'FontName',obj.guiFont,'Units','normalized',...
        'Position',pos,'Min',1,'Max',10000,'FontWeight','bold') ;
    end
    %% table
    function guiCreateTable(obj,name,colnames,colfmt,colwid,coledit,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uitable(parent,'FontName',obj.guiFont,'ColumnName',colnames,'ColumnFormat', colfmt,...
        'ColumnWidth', colwid, 'ColumnEditable',coledit,'Units','normalized','Position',pos);
    end
    %% togglebutton
    function guiCreateTogglebutton(obj,name,txt,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','togglebutton','String',txt,'FontName',obj.guiFont,'Units','normalized',...
        'Position',pos,'FontWeight','bold');
    end
    %% checkbox
    function guiCreateCheckbox(obj,name,txt,val,parent,pos)
      if ischar(parent); parent=obj.gui.(parent); end;
      obj.gui.(name)=uicontrol(parent,'Style','checkbox','String',txt,'FontName',obj.guiFont,...
          'Units','normalized','Position',pos,'FontWeight','bold','Value',val);
    end
  end
  
end

