classdef FlMenu < FlGui
  %FLMENU Application menu for Floodland applications and graphical
  %interface to general Floodland functionality
  %   
  
  properties
    guiTitle='Lucretia:Floodland'; % Title to display on menu GUI 
  end
  properties(Access=private)
    FL % pointer to Floodland object
    appList={}; % list of applications (links to Floodland app objects)
  end
  properties(Dependent)
    appListNames % Names of Floodland apps attached to this menu
  end
  
  methods
    %% constructor
    function obj=FlMenu(FL)
      % Need to pass Floodland object
      if exist('FL','var') && ~strcmp(class(FL),'Floodland')
        error('Must pass Floodland object as first argument');
      end
      obj.FL=FL;
    end
    %% Add application object to menu
    function addApp(obj,appObj)
      % Check requested object to add is a subclass of FlApp
      mc=metaclass(appObj);
      isFlApp=false;
      if ~isempty(mc.SuperclassList)
        for imc=1:length(mc.SuperclassList)
          isFlApp=strcmp(mc.SuperclassList(imc).Name,'FlApp');
          if isFlApp; break; end;
        end
      end
      if ~isFlApp; error('Can only add Floodland application objects to menu, i.e. objects that inherit from FlApp class'); end;
      obj.appList{end+1}=appObj;
    end
    %% Delete application from menu
    function rmApp(obj,appID)
      if exist('appId','var') && isnumeric(appID) && appID>0 && appID<=length(obj.appList)
        obj.appList(appID)=[];
      else
        error('Must supply appID = appList entry (see display(obj))')
      end
    end
    %% Main gui
    function handle=guiMain(obj,~,~)
      if obj.FL.issim
        mode='Sim';
      else
        mode='Live';
      end
      nApp=length(obj.appList);
      border=0.05;
      appButtonSize=50;
      % Generate GUI
      obj.guiCreateFigure('guiMain',sprintf('%s (%s)',obj.guiTitle,mode),[350 max([appButtonSize*nApp,50])]);
      handle=obj.gui.guiMain;
      % Menu & sim mode settings
      mm=uimenu('Parent',obj.gui.guiMain,'Label','Mode');
      obj.gui.menu_simMode=uimenu('Parent',mm,'Label','Sim','Callback',@(src,event)simModeChange(obj,src,event));
      obj.gui.menu_liveMode=uimenu('Parent',mm,'Label','Live','Callback',@(src,event)simModeChange(obj,src,event));
      if obj.FL.issim
        set(obj.gui.menu_simMode,'Checked','on')
        set(obj.gui.menu_liveMode,'Checked','off')
        set(obj.gui.guiMain,'Color',[170,255,128]./255)
      else
        set(obj.gui.menu_simMode,'Checked','off')
        set(obj.gui.menu_liveMode,'Checked','on')
        set(obj.gui.guiMain,'Color',[255,190,0]./255)
      end
%       em=uimenu('Parent',obj.gui.guiMain,'Label','Edit');
      
      % Application buttons
      if ~isempty(obj.appList)
        bh=(1-border*(length(obj.appList)+1))/length(obj.appList);
        for iapp=1:length(obj.appList)
          obj.guiCreatePushbutton(class(obj.appList{iapp}),obj.appList{iapp}.appName,obj.gui.guiMain,...
            [border border+border*(iapp-1)+bh*(iapp-1) 1-2*border bh]);
          set(obj.gui.(class(obj.appList{iapp})),'Callback',@(src,event)guiMain(obj.appList{iapp}))
        end
      end
      drawnow('expose')
    end
    % callback for change mode menu items
    function simModeChange(obj,src,~)
      if src==obj.gui.menu_simMode
        obj.FL.issim=true;
        set(obj.gui.menu_simMode,'Checked','on')
        set(obj.gui.menu_liveMode,'Checked','off')
        set(obj.gui.guiMain,'Color',[170,255,128]./255)
        set(obj.gui.guiMain,'Name',sprintf('%s (%s)',obj.guiTitle,'Sim'))
      else
        obj.FL.issim=false;
        set(obj.gui.menu_simMode,'Checked','off')
        set(obj.gui.menu_liveMode,'Checked','on')
        set(obj.gui.guiMain,'Color',[255,190,0]./255)
        set(obj.gui.guiMain,'Name',sprintf('%s (%s)',obj.guiTitle,'Live'))
      end
      drawnow('expose')
    end
    %% Gets/Sets
    function names=get.appListNames(obj)
      if ~isempty(obj.appList)
        names=cell(1,length(obj.appList));
        for iList=1:length(obj.appList)
          names{iList}=obj.appList{iList}.appName;
        end
      end
    end
    %% Display
    function display(obj)
      % some comment
      names=obj.appListNames;
      fprintf('Index   Application\n')
      fprintf('-----   -----------\n')
      for iapp=1:length(names)
        fprintf('%5d   %s\n',iapp,names{iapp})
      end
    end
  end
  methods(Static)
    %% Status saves
    function statusSave(tag,varargin)
      if ~exist(tag,'dir')
        mkdir(tag);
      end
      % Remove load block if present
      if exist(fullfile(tag,'NOLOAD'),'file')
        delete(fullfile(tag,'NOLOAD'))
      end
      ds=sprintf('%s.mat',datestr(now,30));
      save(fullfile(tag,ds),'varargin')
      
      % Keep max 100 files
      tf=dir(sprintf('%s/*.mat',tag));
      if length(tf)>100
        [a b]=sort({tf.date});
        for ifile=1:length(tf)-100
          delete(fullfile(tag,tf(b(ifile).name)))
        end
      end
    end
    %% Load most recent status save
    function savedata=statusLoad(tag,varargin)
      if ~exist(tag,'dir')
        error('No saved status with this tag')
      end
      % abort if noload tag put on this directory
      if exist(fullfile(tag,'NOLOAD'),'file')
        error('Block put on loading data from this directory, removed with next statusSave')
      end
      tf=dir(sprintf('%s/*.mat',tag));
      [a b]=sort({tf.date});
      sd=load(fullfile(tag,tf(b(end)).name));
      if ~isfield(sd,'varargin') || length(sd.varargin)~=nargin-1
        error('Different number of requested variables than available in saved file')
      end
      for idata=1:length(sd.varargin)
        savedata.(varargin{idata})=sd.varargin{idata};
      end
    end
    %% Block loading from a given tag directory
    function statusLoadBlock(tag)
      evalc(sprintf('!touch %s',fullfile(tag,'NOLOAD')));
    end
  end
end

