classdef FlIndex < handle & FlGui & FlUtils
  %FLINDEX Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Dependent,SetAccess=private)
    PS = []; % PS indexes
    GIRDER = [];
    KLYSTRON = [];
    PS_names % BEAMLINE names of PS's
    PS_indx % BEAMLINE indexes for PS's
    PS_list % Internal PS list elements
    KLYSTRON_names
    KLYSTRON_indx
    KLYSTRON_list
    GIRDER_names
    GIRDER_indx
    GIRDER_list
    NumPars
    KLYSTRON_VMDL % Accelerating potential / MV
    KLYSTRON_PHASEMDL
    PS_BMDL % Ignores PS settings
    INDXnames
    INDXchannels
    INDXbl % BEAMLINE index [upstream middle downstream]
    INDXused % Cell array of names:channels of used controls
  end
  properties(Dependent)
    PS_B % T.m
    PS_KICK % radians
    PS_KMOD % radians / m
    GIRDER_POS
    KLYSTRON_V % Accelerating potential / MV
    KLYSTRON_PHASE % Phase of cavities under this klystron control
    SetPt
    Ampl
    ACT
    Ampl2ACT % conversion ratios between Ampl and ACT units
    limitLow % low control limits
    limitHigh % high control limits
  end
  properties
    useCntrl=[];
    useCntrlChan={};
  end
  properties(Access=protected)
    MasterRef=logical([]);
    MasterInd=[];
    refIndx % Pointer to another FlIndex object instance that has superset of available controls
    refInd % Pointer from this MasterInd to the one in the superset
    indexChoiceFromGui
    indexChanChoiceFromGui
  end
  properties(SetAccess=private)
    ObjID
  end
  properties(Constant)
    Cb = 3.335640952 ; % T.m / GeV
  end
  
  methods
    function obj = FlIndex
      % Constructor
      obj.ObjID=now;
    end
    function addChannel(obj,refInd)
      % Add a new channel to this object from the list available referenced
      % by refIndex
      if refInd<1 || refInd>length(obj.refIndx.MasterInd)
        error('required refInd out of range')
      end
      obj.MasterInd(end+1)=obj.refIndx.MasterInd(refInd);
      obj.MasterRef(end+1,:)=obj.refIndx.MasterRef(refInd,:);
      obj.useCntrl(end+1)=true;
      obj.useCntrlChan{end+1}=obj.refIndx.useCntrlChan{refInd};
      obj.refInd(end+1)=refInd;
      obj.indexChoiceFromGui(end+1)=true;
      obj.indexChanChoiceFromGui{end+1}=obj.useCntrlChan{end};
    end
    function rmChannel(obj,ind)
      % Remove a channel from this object
      if any(ind<1) || any(ind>length(obj.MasterInd))
        error('Required ind out of range')
      end
      obj.MasterInd(ind)=[];
      obj.MasterRef(ind,:)=[];
      obj.useCntrl(ind)=[];
      obj.useCntrlChan(ind)=[];
      obj.refInd(ind)=[];
      obj.indexChoiceFromGui(ind)=[];
      obj.indexChanChoiceFromGui(ind)=[];
    end
  end
  methods
    %% List of used channels
    function chans=get.INDXused(obj)
      irow=0;
      chans={};
      for indx=find(obj.useCntrl)
        if any(obj.useCntrlChan{indx})
          for ichan=find(obj.useCntrlChan{indx})
            irow=irow+1;
            if obj.MasterRef(indx,1)
              channames={'main'};
            elseif obj.MasterRef(indx,2)
              channames={'x' 'dx' 'y' 'dy' 'z' 'dz'};
            elseif obj.MasterRef(indx,3)
              channames={'ampl' 'phase'};
            end
            chans{irow}=sprintf('%s::%s',obj.INDXnames{indx},channames{ichan});
          end
        end
      end
    end
    %% beamline index list
    function indx_out=get.INDXbl(obj)
      global BEAMLINE
      indx=[obj.PS_indx obj.GIRDER_indx obj.KLYSTRON_indx];
      inds=[find(obj.MasterRef(:,1)); find(obj.MasterRef(:,2)); find(obj.MasterRef(:,3))];
      indx=indx(sort(inds));
      indx_out=zeros(length(indx),3);
      for wind=1:length(indx)
        if isfield(BEAMLINE{indx(wind)},'Slices') && ~isempty(BEAMLINE{indx(wind)}.Slices)
          slices=BEAMLINE{indx(wind)}.Slices;
          if length(slices)==1
            indx_out=[indx(wind) indx(wind) indx(wind)];
          else
            indx_out=[slices(1) slices(1)+1 slices(2)];
          end
        else
          indx_out=[indx(wind) indx(wind) indx(wind)];
        end
      end
    end
    %% names list
    function INDXnames=get.INDXnames(obj)
      if ~isempty(obj.MasterInd)
        INDXnames=[obj.PS_names obj.GIRDER_names obj.KLYSTRON_names];
        inds=[find(obj.MasterRef(:,1)); find(obj.MasterRef(:,2)); find(obj.MasterRef(:,3))];
        INDXnames=INDXnames(sort(inds));
      else
        INDXnames={};
      end
    end
    %% define HW
    function defineIndxHW(obj,type,id,pvname,protocol,conv,low,high)
      global PS KLYSTRON GIRDER %#ok<NUSED>
      if ~exist('type','var')
        error('Must supply type (PS, GIRDER or KLYSTRON')
      end
      if ~exist('id','var') || isempty(find(obj.(type)==id, 1))
        error('Must supply type id contained in obj.type')
      end
      if ~exist('pvname','var') || ~isequal(size(pvname),eval(sprintf('size(%s(%d).pvname)',type,id)))
        error('Must supply pvname cell of same dimension as global of same type')
      end
      if ~exist('protocol','var') || ~ismember(protocol,{'AIDA' 'EPICS'})
        error('Must supply protocol as either AIDA or EPICS')
      end
      if ~exist('conv','var') || ~iscell(conv) || ~isequal(size(conv),size(pvname))
        error('Must supply conv conversion scalar or lookup as 2*N array in cell array same dim as pvname')
      end
      if ~exist('low','var') || ~isequal(size(low),size(pvname))
        error('Must supply low bounds for this control (same dimensionality as conv)')
      end
      if ~exist('high','var') || ~isequal(size(high),size(pvname))
        error('Must supply high bounds for this control (same dimensionality as conv)')
      end
      eval(sprintf('%s(%d).pvname=pvname;',type,id));
      eval(sprintf('%s(%d).protocol=protocol;',type,id));
      eval(sprintf('%s(%d).conv=conv;',type,id));
      eval(sprintf('%s(%d).low=low;',type,id));
      eval(sprintf('%s(%d).high=high;',type,id));
      list=obj.(sprintf('%s_list',type));
      ilist=list(obj.(type)==id);
      sz=size(pvname);
      for ipv=1:sz(2)
        obj.useCntrlChan{ilist}(ipv)=~(isempty(pvname{1,ipv}) && isempty(pvname{2,ipv}));
      end
    end
    %% channels list
    function chans=get.INDXchannels(obj)
      global GIRDER KLYSTRON
      if isempty(obj.MasterInd)
        chans={};
        return
      end
      chans=cell(2,length(obj.MasterInd));
      for rw=1:2
        for ind=1:length(obj.MasterInd)
          if obj.MasterRef(ind,1)
            chans{rw,ind}=1;
          elseif obj.MasterRef(ind,2)
            for ichan=1:6
              if ~isempty(GIRDER{obj.MasterInd(ind)}.pvname{rw,ichan})
                chans{rw,ind}(length(chans{rw,ind})+1)=ichan;
              end
            end
          elseif obj.MasterRef(ind,3)
            for ichan=1:2
              if ~isempty(KLYSTRON(obj.MasterInd(ind)).pvname{rw,ichan})
                chans{rw,ind}(length(chans{rw,ind})+1)=ichan;
              end
            end
          end
        end
      end
    end
    %% high/low limits
    function vals=get.limitLow(obj)
      global PS GIRDER KLYSTRON
      pslist=obj.PS_list;
      psind=obj.PS;
      girlist=obj.GIRDER_list;
      girind=obj.GIRDER;
      klylist=obj.KLYSTRON_list;
      klyind=obj.KLYSTRON;
      vals=cell(1,length(obj.MasterInd));
      if ~isempty(pslist)
        for ips=1:length(pslist)
          vals{pslist(ips)}=PS(psind(ips)).low;
        end
      end
      if ~isempty(girlist)
        for igir=1:length(girlist)
          vals{girlist(igir)}=GIRDER{girind(igir)}.low;
        end
      end
      if ~isempty(klylist)
        for ikly=1:length(klylist)
          vals{klylist(ikly)}=KLYSTRON(klyind(ikly)).low;
        end
      end
    end
    function vals=get.limitHigh(obj)
      global PS GIRDER KLYSTRON
      pslist=obj.PS_list;
      psind=obj.PS;
      girlist=obj.GIRDER_list;
      girind=obj.GIRDER;
      klylist=obj.KLYSTRON_list;
      klyind=obj.KLYSTRON;
      vals=cell(1,length(obj.MasterInd));
      if ~isempty(pslist)
        for ips=1:length(pslist)
          vals{pslist(ips)}=PS(psind(ips)).high;
        end
      end
      if ~isempty(girlist)
        for igir=1:length(girlist)
          vals{girlist(igir)}=GIRDER{girind(igir)}.high;
        end
      end
      if ~isempty(klylist)
        for ikly=1:length(klylist)
          vals{klylist(ikly)}=KLYSTRON(klyind(ikly)).high;
        end
      end
    end
    function vals=get.Ampl2ACT(obj)
      % Scale factor relating obj.Ampl values to obj.ACT values
      % obj.ACT{icntrl}(ichan)=obj.Ampl{icntrl}(ichan)*vals{icntrl}(ichan)
      curAmpl=obj.Ampl;
      newAmpl=curAmpl;
      vals=newAmpl;
      for icntrl=1:length(curAmpl)
        for ichan=1:length(curAmpl{icntrl})
          newAmpl{icntrl}(ichan)=1;
        end
      end
      obj.Ampl=newAmpl;
      act=obj.ACT;
      for icntrl=1:length(curAmpl)
        for ichan=1:length(curAmpl{icntrl})
          vals{icntrl}(ichan)=act{icntrl}(ichan);
        end
      end
      obj.Ampl=curAmpl;
    end
    %% SetPt list
    function vals=get.SetPt(obj)
      global PS GIRDER KLYSTRON
      pslist=obj.PS_list;
      psind=obj.PS;
      girlist=obj.GIRDER_list;
      girind=obj.GIRDER;
      klylist=obj.KLYSTRON_list;
      klyind=obj.KLYSTRON;
      vals=cell(1,length(obj.MasterInd));
      if ~isempty(pslist)
        for ips=1:length(pslist)
          vals{pslist(ips)}=PS(psind(ips)).SetPt;
        end
      end
      if ~isempty(girlist)
        for igir=1:length(girlist)
          vals{girlist(igir)}=GIRDER{girind(igir)}.MoverSetPt;
        end
      end
      if ~isempty(klylist)
        for ikly=1:length(klylist)
          vals{klylist(ikly)}=[KLYSTRON(klyind(ikly)).AmplSetPt KLYSTRON(klyind(ikly)).PhaseSetPt];
        end
      end
    end
    function set.SetPt(obj,vals)
      global PS GIRDER KLYSTRON
      pslist=obj.PS_list;
      psind=obj.PS;
      girlist=obj.GIRDER_list;
      girind=obj.GIRDER;
      klylist=obj.KLYSTRON_list;
      klyind=obj.KLYSTRON;
      if ~isempty(pslist)
        for ips=1:length(pslist)
          PS(psind(ips)).SetPt=vals{pslist(ips)};
        end
      end
      if ~isempty(girlist)
        for igir=1:length(girlist)
          GIRDER{girind(igir)}.MoverSetPt=vals{girlist(igir)};
        end
      end
      if ~isempty(klylist)
        for ikly=1:length(klylist)
          KLYSTRON(klyind(ikly)).AmplSetPt=vals{klylist(ikly)}(1);
          KLYSTRON(klyind(ikly)).PhaseSetPt=vals{klylist(ikly)}(2);
        end
      end
    end
    function vals=get.Ampl(obj)
      % Get Ampl readback values (PS.Ampl / GIRDER.MoverSetPt /
      % KLYSTRON.Ampl / KLYSTRON.Phase)
      global PS GIRDER KLYSTRON
      pslist=obj.PS_list;
      psind=obj.PS;
      girlist=obj.GIRDER_list;
      girind=obj.GIRDER;
      klylist=obj.KLYSTRON_list;
      klyind=obj.KLYSTRON;
      vals=cell(1,length(obj.MasterInd));
      if ~isempty(pslist)
        for ips=1:length(pslist)
          vals{pslist(ips)}=PS(psind(ips)).Ampl;
        end
      end
      if ~isempty(girlist)
        for igir=1:length(girlist)
          vals{girlist(igir)}=GIRDER{girind(igir)}.MoverPos;
        end
      end
      if ~isempty(klylist)
        for ikly=1:length(klylist)
          vals{klylist(ikly)}=[KLYSTRON(klyind(ikly)).Ampl KLYSTRON(klyind(ikly)).Phase];
        end
      end
    end
    function set.Ampl(obj,vals)
      % Set Ampl readback values (PS.Ampl / GIRDER.MoverSetPt /
      % KLYSTRON.Ampl / KLYSTRON.Phase) Note that this will set the
      % READBACK values, and this will be overwritten by a control system
      % get command (Floodland hwGet method).
      global PS GIRDER KLYSTRON
      pslist=obj.PS_list;
      psind=obj.PS;
      girlist=obj.GIRDER_list;
      girind=obj.GIRDER;
      klylist=obj.KLYSTRON_list;
      klyind=obj.KLYSTRON;
      if ~isempty(pslist)
        for ips=1:length(pslist)
          PS(psind(ips)).Ampl=vals{pslist(ips)};
        end
      end
      if ~isempty(girlist)
        for igir=1:length(girlist)
          GIRDER{girind(igir)}.MoverPos=vals{girlist(igir)};
        end
      end
      if ~isempty(klylist)
        for ikly=1:length(klylist)
          KLYSTRON(klyind(ikly)).Ampl=vals{klylist(ikly)}(1);
          KLYSTRON(klyind(ikly)).Phase=vals{klylist(ikly)}(2);
        end
      end
    end
    function vals=get.ACT(obj)
      % Get controls values in Physical Lucretia units (e.g. T.m for Quads
      % and Volts for klystron Amplutudes)
      pslist=obj.PS_list;
      girlist=obj.GIRDER_list;
      klylist=obj.KLYSTRON_list;
      vals=cell(1,length(obj.MasterInd));
      if ~isempty(pslist)
        psb=obj.PS_B;
        for ips=1:length(pslist)
          vals{pslist(ips)}=psb(ips);
        end
      end
      if ~isempty(girlist)
        gpos=obj.GIRDER_POS;
        for igir=1:length(girlist)
          vals{girlist(igir)}=gpos(igir,:);
        end
      end
      if ~isempty(klylist)
        kv=obj.KLYSTRON_V;
        kp=obj.KLYSTRON_PHASE;
        for ikly=1:length(klylist)
          vals{klylist(ikly)}=[kv(ikly) kp(ikly)];
        end
      end
    end
    function set.ACT(obj,vals)
      % Set controls values in Physical Lucretia units (e.g. T.m for Quads
      % and Volts for klystron Amplutudes)
      pslist=obj.PS_list;
      girlist=obj.GIRDER_list;
      klylist=obj.KLYSTRON_list;
      if ~isempty(pslist)
        psb=obj.PS_B;
        for ips=1:length(pslist)
          psb(ips)=vals{pslist(ips)};
        end
        obj.PS_B=psb;
      end
      if ~isempty(girlist)
        gpos=obj.GIRDER_POS;
        for igir=1:length(girlist)
          gpos(igir,:)=vals{girlist(igir)};
        end
        obj.GIRDER_POS=gpos;
      end
      if ~isempty(klylist)
        kv=obj.KLYSTRON_V;
        kp=obj.KLYSTRON_PHASE;
        for ikly=1:length(klylist)
          kv(ikly)=vals{klylist(ikly)}(1);
          kp(ikly)=vals{klylist(ikly)}(2);
        end
        obj.KLYSTRON_V=kv;
        obj.KLYSTRON_PHASE=kp;
      end
    end
    %% PS list
    function list=get.PS_list(obj)
      list=find(obj.MasterRef(:,1));
    end
     %% GIRDER list
    function list=get.GIRDER_list(obj)
      list=find(obj.MasterRef(:,2));
    end
     %% KLYSTRON list
    function list=get.KLYSTRON_list(obj)
      list=find(obj.MasterRef(:,3));
    end
    %% GIRDER POS get
    function gpos=get.GIRDER_POS(obj)
      global GIRDER
      gpos=cell2mat(arrayfun(@(x) GIRDER{x}.MoverPos,obj.GIRDER,'UniformOutput',false));
    end
    %% GIRDER POS set
    function set.GIRDER_POS(obj,value)
      global GIRDER
      for igir=obj.GIRDER
        GIRDER{igir}.MoverSetPt=value(igir,:);
      end
    end
    %% PS list
    function pslist=get.PS(obj)
      pslist=obj.MasterInd(obj.MasterRef(:,1));
    end
    %% GIRDER list
    function glist=get.GIRDER(obj)
      glist=obj.MasterInd(obj.MasterRef(:,2));
    end
    %% KLYSTRON list
    function klist=get.KLYSTRON(obj)
      klist=obj.MasterInd(obj.MasterRef(:,3));
    end
    %% Get KLYSTRON V
    function v=get.KLYSTRON_V(obj)
      global KLYSTRON BEAMLINE
      for ikly=obj.KLYSTRON
        v(ikly)=KLYSTRON(obj.KLYSTRON(ikly)).Ampl*sum(arrayfun(@(x) BEAMLINE{x}.Volt,KLYSTRON(obj.KLYSTRON(ikly)).Element));
      end
    end
    %% Set KLYSTRON V
    function set.KLYSTRON_V(obj,value)
      global KLYSTRON BEAMLINE
      for ikly=obj.KLYSTRON
        KLYSTRON(obj.KLYSTRON(ikly)).AmplSetPt=value(ikly)/sum(arrayfun(@(x) BEAMLINE{x}.Volt,KLYSTRON(obj.KLYSTRON(ikly)).Element));
      end
    end
    %% Get KLYSTRON PHASE
    function ph=get.KLYSTRON_PHASE(obj)
      global KLYSTRON BEAMLINE
      for ikly=obj.KLYSTRON
        ph(ikly)=KLYSTRON(obj.KLYSTRON(ikly)).Phase+BEAMLINE{KLYSTRON(obj.KLYSTRON(ikly)).Element(1)}.Phase;
      end
    end
    %% Set KLYSTRON PHASE
    function set.KLYSTRON_PHASE(obj,value)
      global KLYSTRON BEAMLINE
      for ikly=obj.KLYSTRON
        KLYSTRON(obj.KLYSTRON(ikly)).PhaseSetPt=value(ikly)-BEAMLINE{KLYSTRON(obj.KLYSTRON(ikly)).Element(1)}.Phase;
      end
    end
     %% Get KLYSTRON V MDL
    function v=get.KLYSTRON_VMDL(obj)
      global KLYSTRON BEAMLINE
      for ikly=obj.KLYSTRON
        v(ikly)=sum(arrayfun(@(x) BEAMLINE{x}.Volt,KLYSTRON(obj.KLYSTRON(ikly)).Element));
      end
    end
    %% Get KLYSTRON PHASE MDL
    function ph=get.KLYSTRON_PHASEMDL(obj)
      global KLYSTRON BEAMLINE
      for ikly=obj.KLYSTRON
        ph(ikly)=BEAMLINE{KLYSTRON(obj.KLYSTRON(ikly)).Element(1)}.Phase;
      end
    end
    %% Get NumPars
    function numPars=get.NumPars(obj)
      numPars=length(obj.PS)+length(obj.KLYSTRON)+length(obj.GIRDER);
    end
    %% Set corrector KMOD
    function set.PS_KMOD(obj,val)
      global BEAMLINE
      obj.PS_KICK=val.*arrayfun(@(x) BEAMLINE{x}.L(1),obj.PS_indx);
    end
    %% Get corrector KMOD
    function ret = get.PS_KMOD(obj)
      global BEAMLINE
      ret=obj.PS_KICK./arrayfun(@(x) BEAMLINE{x}.L(1),obj.PS_indx);
    end
    %% Set corrector KICK
    function set.PS_KICK(obj,val)
      global BEAMLINE
      bind=obj.PS_indx;
      corb=zeros(size(obj.PS));
      for icor=1:length(bind)
        corb(icor)=obj.Cb * BEAMLINE{bind(icor)}.P.*val(icor);
      end
      obj.PS_B=corb;
    end
    %% Get corrector KICK
    function ret = get.PS_KICK(obj)
      global BEAMLINE
      bind=obj.PS_indx;
      corb=obj.PS_B;
      kick=zeros(size(obj.PS));
      for icor=1:length(bind)
        kick(icor)=corb(icor)./(obj.Cb * BEAMLINE{bind(icor)}.P);
      end
      ret=kick;
    end
    %% Get corrector B
    function ret = get.PS_B(obj)
      global BEAMLINE PS
      ret=arrayfun(@(x) PS(x).Ampl.*sum(arrayfun(@(xx) BEAMLINE{xx}.B,PS(x).Element)),obj.PS);
    end
     %% Get corrector B (Model)
    function ret = get.PS_BMDL(obj)
      global BEAMLINE PS
      ret=arrayfun(@(x) sum(arrayfun(@(xx) BEAMLINE{xx}.B,PS(x).Element)),obj.PS);
    end
    %% Set corrector B
    function set.PS_B(obj,val)
      global BEAMLINE PS
      for ips=obj.PS
        PS(ips).SetPt=val(ips)./sum(arrayfun(@(x) BEAMLINE{x}.B,PS(ips).Element));
      end
    end
    %% Get Klystron name list
    function names=get.KLYSTRON_names(obj)
      global BEAMLINE
      names=arrayfun(@(x) BEAMLINE{x}.Name,obj.KLYSTRON_indx,'UniformOutput',false);
    end
    %% Get GIRDER name list
    function names=get.GIRDER_names(obj)
      global BEAMLINE
      names=arrayfun(@(x) BEAMLINE{x}.Name,obj.GIRDER_indx,'UniformOutput',false);
    end
    %% Get PS name list
    function names=get.PS_names(obj)
      global BEAMLINE
      names=arrayfun(@(x) BEAMLINE{x}.Name,obj.PS_indx,'UniformOutput',false);
    end
    %% Get list of PS by BEAMLINE indices
    function psind=get.PS_indx(obj)
      global PS
      psind=arrayfun(@(x) PS(x).Element(1),obj.PS);
    end
    %% Get list of KLYS by BEAMLINE indices
    function klyind=get.KLYSTRON_indx(obj)
      global KLYSTRON
      klyind=arrayfun(@(x) KLYSTRON(x).Element(1),obj.KLYSTRON);
    end
    %% Get list of GIRDER by BEAMLINE indices
    function girind=get.GIRDER_indx(obj)
      global GIRDER
      girind=arrayfun(@(x) GIRDER{x}.Element(1),obj.GIRDER);
    end
    %% Display
    function display(obj)
      global BEAMLINE KLYSTRON PS GIRDER %#ok<NUSED>
       fprintf('====================================\n')
       fprintf('   ID        Type          Name\n')
       fprintf('====================================\n')
       t={'PS' 'GIRDER' 'KLYSTRON'};
       if ~isempty(obj.MasterRef)
         for id=1:length(obj.MasterRef(:,1))
           type=t{obj.MasterRef(id,:)};
           str=sprintf('%s(%d)',type,obj.MasterInd(id));
           if isempty(BEAMLINE)
             name='??';
           else
             try
               evalc(sprintf('name=BEAMLINE{%s.Element(1)}.Name',str));
             catch
               name='??';
             end
           end
           fprintf('%5d   %15s    %s\n',id,str,name)
         end
       end
    end
    %% +
    function plus(obj,A)
      mc=metaclass(A);
      isFlIndex=strcmp(class(A),'FlIndex');
      if ~isempty(mc.SuperclassList) && ~isFlIndex
        for imc=1:length(mc.SuperclassList)
          isFlIndex=strcmp(mc.SuperclassList(imc).Name,'FlIndex');
          if isFlIndex; break; end;
        end
      end
      if ~isFlIndex; error('Addition only supported for 2 classes that are or inherit from FlIndex'); end;
      aindref=find(ismember(A.MasterInd,obj.MasterInd));
      aref=true(size(A.MasterInd));
      for iref=1:length(aindref)
        oref=find(obj.MasterInd==A.MasterInd(aindref(iref)));
        for ioref=1:length(oref)
          if isequal(A.MasterRef(aindref(iref),:),obj.MasterRef(oref(ioref),:))
            aref(aindref(iref))=false;
          end
        end
      end
      obj.MasterInd=[obj.MasterInd; A.MasterInd(aref)];
      obj.MasterRef=[obj.MasterRef; A.MasterRef(aref,:)];
      obj.useCntrl=[obj.useCntrl A.useCntrl(aref)];
      obj.useCntrlChan(end+1:end+sum(aref))=A.useCntrlChan(aref);
    end
    %% -
    function minus(obj,A)
      mc=metaclass(A);
      isFlIndex=strcmp(class(A),'FlIndex');
      if ~isempty(mc.SuperclassList) && ~isFlIndex
        for imc=1:length(mc.SuperclassList)
          isFlIndex=strcmp(mc.SuperclassList(imc).Name,'FlIndex');
          if isFlIndex; break; end;
        end
      end
      if ~isFlIndex; error('Addition only supported for 2 classes that are or inherit from FlIndex'); end;
      aindref=find(ismember(A.MasterInd,obj.MasterInd));
      for iref=1:length(aindref)
        oref=find(obj.MasterInd==A.MasterInd(aindref(iref)));
        for ioref=1:length(oref)
          if isequal(A.MasterRef(aindref(iref),:),obj.MasterRef(oref(ioref),:))
            obj.MasterInd(oref(ioref))=[];
            obj.MasterRef(oref(ioref),:)=[];
          end
        end
      end
    end
    %% ( ... )
%     function subsref(obj,S)
%       switch S(1).type
%         case '.'
%           builtin('subsref',obj,S)
%         case '()'
%           B=FlIndex;
%           if length(S)>1 && strcmp(S(2).type,'.') && ismember(S(2).subs,{'PS','GIRDER','KLYSTRON'})
%             t={'PS' 'GIRDER' 'KLYSTRON'};
%             tind=ismember(obj.MasterInd(obj.MasterRef(:,ismember(t,S(2).subs{1}))),S(1).subs{1});
%             B.MasterRef=obj.MasterRef(tind,:);
%             B.MasterInd=obj.MasterInd(tind,:);
%           else
%             B.MasterRef=obj.MasterRef(S(1).subs{1},:);
%             B.MasterInd=obj.MasterInd(S(1).subs{1});
%           end
%         case '{}'
%           error('FlIndex:subsref','Not a supported subscripted reference')
%       end
%     end
    %% Add PS
    function obj = addPS(obj, blList)
      global BEAMLINE PS
      if ~iscell(blList); blList={blList}; end;
      % If no PS exist for this index, first make it
      for ibl=1:length(blList)
        if ~isfield(BEAMLINE{blList{ibl}(1)},'PS') || ~BEAMLINE{blList{ibl}(1)}.PS
          stat = AssignToPS( blList{ibl}, length(PS)+1 ) ;
          if stat{1}~=1; error('PS assignment error:\n%s\n',stat{2}); end;
          ind=length(PS);
        else
          ind=BEAMLINE{blList{ibl}(1)}.PS;
          if ismember(ind,obj.PS); continue; end;
        end
        if ~isfield(PS(ind),'pvname') || ~iscell(PS(ind).pvname)
          PS(ind).trimpv=cell(2,1);
          PS(ind).trimconv={1;1};
          PS(ind).protocol='none';
          PS(ind).pvname=cell(2,1);
          PS(ind).preCommand=cell(2,1);
          PS(ind).postCommand=cell(2,1);
          PS(ind).conv={1;1};
          PS(ind).unipolar=false;
          PS(ind).timestamp=0;
          PS(ind).high=NaN(2,1);
          PS(ind).low=NaN(2,1);
          PS(ind).trimHigh=NaN(2,1);
          PS(ind).trimLow=NaN(2,1);
          PS(ind).trimUnipolar=false;
          obj.MasterInd(end+1)=ind;
          obj.MasterRef(end+1,:)=[true false false];
          obj.useCntrl(end+1)=true;
          obj.useCntrlChan{end+1}=true;
          obj.indexChoiceFromGui(end+1)=true;
          obj.indexChanChoiceFromGui{end+1}=obj.useCntrlChan{end};
        end
      end
    end
    %% Add Klystrons
    function obj = addKlystron(obj, kList)
      global BEAMLINE KLYSTRON
      if ~iscell(kList); kList={kList}; end;
      % If no KLYSTRON exist for this index, first make it
      for ik=1:length(kList)
        if ~isfield(BEAMLINE{kList{ik}(1)},'Klystron') || ~BEAMLINE{kList{ik}(1)}.Klystron
          stat = AssignToKlystron(kList{ik}, length(KLYSTRON)+1) ;
          if stat{1}~=1; error('KLYSTRON assignment error:\n%s\n',stat{2}); end;
          ind=length(KLYSTRON);
        else
          ind=BEAMLINE{kList{ik}(1)}.Klystron;
          if ismember(ind,obj.KLYSTRON); continue; end;
        end
        if ~isfield(KLYSTRON(ind),'pvname') || ~iscell(KLYSTRON(ind).pvname)
          KLYSTRON(ind).cunits='none';
          KLYSTRON(ind).protocol='none';
          KLYSTRON(ind).pvname=cell(2,2);
          KLYSTRON(ind).preCommand=cell(2,2);
          KLYSTRON(ind).postCommand=cell(2,2);
          KLYSTRON(ind).conv={1 1;1 1};
          KLYSTRON(ind).timestamp=0;
          KLYSTRON(ind).high=NaN(2,2);
          KLYSTRON(ind).low=NaN(2,2);
          KLYSTRON(ind).statpv=[];
        end
        obj.MasterInd(end+1)=ind;
        obj.MasterRef(end+1,:)=[false false true];
        obj.useCntrl(end+1)=true;
        obj.useCntrlChan{end+1}=[true true];
        obj.indexChoiceFromGui(end+1)=true;
        obj.indexChanChoiceFromGui{end+1}=obj.useCntrlChan{end};
      end
    end
    %% Add Movers
    function obj = addMover(obj, mList)
      global BEAMLINE GIRDER
      % If no GIRDER exist for this index, first make it
      for im=1:length(mList)
        if ~isfield(BEAMLINE{mList(im)},'Girder') || ~BEAMLINE{mList(im)}.Girder
          stat = AssignToGirder(mList(im), length(GIRDER)+1, 1) ;
          if stat{1}~=1; error('GIRDER assignment error:\n%s\n',stat{2}); end;
          ind=length(GIRDER);
        else
          ind=BEAMLINE{mList(im)}.Girder;
        end
        if ~isfield(GIRDER{ind},'pvname') || ~iscell(GIRDER{ind}.pvname)
          GIRDER{ind}.Mover=1:6;
          GIRDER{ind}.cunits=cell(2,6);
          GIRDER{ind}.protocol=cell(2,6);
          GIRDER{ind}.pvname=cell(2,6);
          GIRDER{ind}.preCommand=cell(2,6);
          GIRDER{ind}.postCommand=cell(2,6);
          GIRDER{ind}.conv={1 1 1 1 1 1; 1 1 1 1 1 1};
          GIRDER{ind}.timestamp=0;
          GIRDER{ind}.high=NaN(2,6);
          GIRDER{ind}.low=NaN(2,6);
        end
        obj.MasterInd(end+1)=ind;
        obj.MasterRef(end+1,:)=[false true false];
        obj.useCntrl(end+1)=true;
        obj.useCntrlChan{end+1}=[true true true true true true];
        obj.indexChoiceFromGui(end+1)=true;
        obj.indexChanChoiceFromGui{end+1}=obj.useCntrlChan{end};
      end
    end
  end
  %% GUI
  methods
    % Process GUI callbacks
    function guiIndexCallback(obj,src,~)
      try
        if src==obj.gui.display_s
          if get(obj.gui.display_s,'Value')
            set(obj.gui.display_alpha,'Value',false)
            obj.updateIndexGui;
          else
            set(obj.gui.display_s,'Value',true)
          end
        elseif src==obj.gui.display_alpha
          if get(obj.gui.display_alpha,'Value')
            set(obj.gui.display_s,'Value',false)
            obj.updateIndexGui;
          else
            set(obj.gui.display_alpha,'Value',true)
          end
        elseif src==obj.gui.select_accept
          delete(obj.gui.guiIndexChoice);
          obj.gui.guiIndexChoice=[];
        elseif src==obj.gui.select_cancel || src==obj.gui.guiIndexChoice
          obj.indexChoiceFromGui=[];
          obj.indexChanChoiceFromGui=[];
          delete(obj.gui.guiIndexChoice);
          obj.gui.guiIndexChoice=[];
        elseif src==obj.gui.selbutton1
          sel=get(obj.gui.indexCbox1,'Value');
          mI=get(obj.gui.indexCbox1,'UserData');
          obj.indexChoiceFromGui(mI(sel))=false;
          obj.indexChanChoiceFromGui(mI(sel))=obj.useCntrlChan(mI(sel));
          set(obj.gui.indexCbox1,'Value',1)
          obj.updateIndexGui;
        elseif src==obj.gui.selbutton2
          sel=get(obj.gui.indexCbox2,'Value');
          mI=get(obj.gui.indexCbox2,'UserData');
          obj.indexChoiceFromGui(mI(sel))=true;
          obj.indexChanChoiceFromGui(mI(sel))=obj.useCntrlChan(mI(sel));
          set(obj.gui.indexCbox2,'Value',1)
          obj.updateIndexGui;
        else
          obj.updateIndexGui;
        end
      catch ME
        disp(ME.message)
        delete(gcf)
      end
    end
    % Update GUI fields
    function updateIndexGui(obj)
      global PS GIRDER KLYSTRON BEAMLINE %#ok<NUSED>
      % Get display options
      cname={'PS' 'GIRDER' 'KLYSTRON'};
      showControls=[isfield(obj.gui,'controls_PS')&&get(obj.gui.controls_PS,'Value') ...
                    isfield(obj.gui,'controls_GIRDER')&&get(obj.gui.controls_GIRDER,'Value') ...
                    isfield(obj.gui,'controls_KLYSTRON')&&get(obj.gui.controls_KLYSTRON,'Value')];
      if ~any(showControls); return; end;
      alltypes={'QUAD' 'SEXT' 'MULT' 'COR'};
      types=alltypes(logical([get(obj.gui.type_QUAD,'Value') get(obj.gui.type_SEXT,'Value') get(obj.gui.type_MULT,'Value') get(obj.gui.type_COR,'Value')]));
      if ismember('COR',types); types{end}='XCOR'; types{end+1}='YCOR'; end;
      % Form right display
      displayStr2={}; s2=[]; name2={}; mI=[];
      chnames={{'main'} {'x' 'dx' 'y' 'dy' 'z' 'dz'} {'ampl' 'pha'}};
      allchans={};
      chansWithHW=obj.INDXchannels;
      for ic=find(showControls)
        if isfield(obj.gui,sprintf('controls_%s',cname{ic})) && ~isempty(obj.(cname{ic}))
          blInd=obj.([cname{ic} '_indx']);
          typeInd=obj.(cname{ic});
          masterInd=obj.([cname{ic} '_list']);
          for ibl=1:length(blInd)
            if ismember(BEAMLINE{blInd(ibl)}.Class,types) || ~isempty(regexp(BEAMLINE{blInd(ibl)}.Class,'CAV$', 'once'))
              thisChanStr='(';
              chans=find(obj.indexChanChoiceFromGui{masterInd(ibl)});
              if ~isempty(chans)
                for ichan=chans
                  if ismember(ichan,chansWithHW{1,masterInd(ibl)}) || ismember(ichan,chansWithHW{2,masterInd(ibl)})
                    thisChanStr=[thisChanStr chnames{ic}{ichan} ','];
                    allchans{end+1}=chnames{ic}{ichan};
                  end
                end
              end
              thisChanStr=[thisChanStr ')']; thisChanStr=regexprep(thisChanStr,',)',')');
              displayStr2{end+1}=sprintf('%d: %s [%s(%d)] %s',blInd(ibl),BEAMLINE{blInd(ibl)}.Name,cname{ic},typeInd(ibl),thisChanStr) ;
              s2(end+1)=blInd(ibl);
              name2{end+1}=BEAMLINE{blInd(ibl)}.Name;
              mI(end+1)=masterInd(ibl);
            end
          end
        end
      end
      sel=get(obj.gui.indexCbox2,'Value');
      % Set channel selection list with options
      set(obj.gui.chansel,'Value',1);
      set(obj.gui.chansel,'String',unique(allchans));
      % - put in required sort order
      if max(sel)>length(displayStr2)
        set(obj.gui.indexCbox2,'Value',1)
      end
      name1=name2(ismember(mI,find(obj.indexChoiceFromGui)));
      s1=s2(ismember(mI,find(obj.indexChoiceFromGui)));
      name2=name2(~ismember(mI,find(obj.indexChoiceFromGui)));
      s2=s2(~ismember(mI,find(obj.indexChoiceFromGui)));
      displayStr1=displayStr2(ismember(mI,find(obj.indexChoiceFromGui)));
      mI1=mI(ismember(mI,find(obj.indexChoiceFromGui)));
      displayStr2=displayStr2(~ismember(mI,find(obj.indexChoiceFromGui)));
      mI2=mI(~ismember(mI,find(obj.indexChoiceFromGui)));
      if get(obj.gui.display_alpha,'Value')
        [Y I]=sort(name2);
        [Y I1]=sort(name1);
      else
        [Y I]=sort(s2);
        [Y I1]=sort(s1);
      end
      set(obj.gui.indexCbox2,'String',displayStr2(I))
      set(obj.gui.indexCbox2,'UserData',mI2(I))
      set(obj.gui.indexCbox1,'String',displayStr1(I1))
      set(obj.gui.indexCbox1,'UserData',mI1(I1))
      drawnow('expose');
    end
    % Setup GUI and display
    function han=guiIndexChoice(obj)
      % If nothing defined just exit
      if isempty(obj.MasterInd)
        han=[];
        return
      end
      % Choose categories of controls to choose from
      hwChoice={'PS' 'GIRDER' 'KLYSTRON'};
      % Create main figure window
      obj.guiCreateFigure('guiIndexChoice','Choose Controls',[800 600]);
      set(obj.gui.guiIndexChoice,'CloseRequestFcn',@(src,event)guiIndexCallback(obj,src,event));
      han=obj.gui.guiIndexChoice;
      % Define borders
      border_top=0.02; 
      border_bottom=0.02;
      border_left=0.02;
      border_right=0.02;
      % Define element areas
      area_ver=[80 5 15].*0.01; %
      area_hor={[45 10 45].*0.01 [22 4 22 4 22 4 22].*0.01};
      % User Area
      userSize=[(1-(border_left+border_right)) (1-(border_top+border_bottom))];
      % Define choice box panels
      cpanpos{1}=[border_left border_bottom+(area_ver(3)+area_ver(2))*userSize(2)];
      cpanpos{2}=[border_left+(area_hor{1}(1)+area_hor{1}(2))*userSize(1) ...
        border_bottom+(area_ver(3)+area_ver(2))*userSize(2)];
      cpansize(1)=userSize(1)*area_hor{1}(1) ;
      cpansize(2)=userSize(2)*area_ver(1) ;
      cpanborder=[0.01 0.01];
      obj.guiCreatePanel('cbox1_panel','Controls Selection','guiIndexChoice',[cpanpos{1}(1) cpanpos{1}(2) cpansize]);
      obj.guiCreatePanel('cbox2_panel','Controls Available','guiIndexChoice',[cpanpos{2}(1) cpanpos{2}(2) cpansize]);
      % Define choice listboxes
      obj.guiCreateListbox('indexCbox1','','cbox1_panel',[cpanborder 1-2*cpanborder]);
      set(obj.gui.indexCbox1,'Callback',@(src,event)guiCselCallback(obj,src,event));
      obj.guiCreateListbox('indexCbox2','','cbox2_panel',[cpanborder 1-2*cpanborder]);
      % Define selection buttons
      selborder=0.01;
      selsize=[area_hor{1}(2)*userSize(1)*0.9 cpansize(2)*0.1];
      selpos=[border_left+userSize(1)*(area_hor{1}(1)+area_hor{1}(2)*selborder) ...
        border_bottom+userSize(2)*(area_ver(3)+area_ver(2))]+[0 cpansize(2)*0.6];
      obj.guiCreatePushbutton('selbutton1','->','guiIndexChoice',[selpos selsize]);
      set(obj.gui.selbutton1,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      selpos=selpos-[0 cpansize(2)*0.2];
      obj.guiCreatePushbutton('selbutton2','<-','guiIndexChoice',[selpos selsize]);
      set(obj.gui.selbutton2,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      % Define option panels
      oppansize{1}=[area_hor{2}(1).*userSize(1) area_ver(3)*userSize(2)];
      oppansize{2}=[area_hor{2}(3).*userSize(1) area_ver(3)*userSize(2)];
      oppansize{3}=[area_hor{2}(5).*userSize(1) area_ver(3)*userSize(2)];
      oppansize{4}=[area_hor{2}(7).*userSize(1) area_ver(3)*userSize(2)];
      osize=(1-6*border_left)/5;
      obj.guiCreatePanel('oppan1','Sort by','guiIndexChoice',[border_left border_bottom osize oppansize{1}(2)]);
      obj.guiCreatePanel('oppan1b','Channels','guiIndexChoice',[border_left*2+osize border_bottom osize oppansize{1}(2)]);
      obj.guiCreatePanel('oppan2','Show Controls','guiIndexChoice',[border_left*3+osize*2 border_bottom osize oppansize{1}(2)]);
      obj.guiCreatePanel('oppan3','Type','guiIndexChoice',[border_left*4+osize*3 border_bottom osize oppansize{1}(2)]);
      obj.guiCreatePanel('oppan4','Make Selection','guiIndexChoice',[border_left*5+osize*4 border_bottom osize oppansize{1}(2)]);
      % List options
      % - display
      border=0.02;
      osize=(1-3*border)/2;
      obj.guiCreateRadiobutton('display_alpha','Name',0,'oppan1',[border border 1-2*border osize])
      set(obj.gui.display_alpha,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      obj.guiCreateRadiobutton('display_s','S',1,'oppan1',[border border*2+osize 1-2*border osize])
      set(obj.gui.display_s,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      % - Channels
      border=0.02;
      osize=(1-3*border)/2;
      obj.guiCreateTogglebutton('showchan','Select Channel(s)','oppan1b',[border border 1-2*border osize]);
      set(obj.gui.showchan,'Callback',@(src,event)guiIndexChanSelCallback(obj,src,event));
      obj.guiCreatePopupmenu('chansel','Toggle selection:','oppan1b',[border border*2+osize 1-2*border osize]);
      set(obj.gui.chansel,'Callback',@(src,event)guiChanSelCallback(obj,src,event));
      set(obj.gui.chansel,'Value',1);
      % - show controls
      nchoice=length(hwChoice);
      osize=(1-(1+nchoice)*border)/nchoice;
      for ic=1:nchoice
        obj.guiCreateCheckbox(sprintf('controls_%s',hwChoice{ic}),hwChoice{ic},1,'oppan2',[border border*ic+osize*(ic-1) 1-2*border osize]);
        set(obj.gui.(sprintf('controls_%s',hwChoice{ic})),'Callback',@(src,event)guiIndexCallback(obj,src,event)) ;
      end
      % - type selection
      osize=(1-5*border)/4;
      obj.guiCreateCheckbox('type_QUAD','QUAD',1,'oppan3',[border border 1-2*border osize]);
      set(obj.gui.type_QUAD,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      obj.guiCreateCheckbox('type_SEXT','SEXT',1,'oppan3',[border border*2+osize 1-2*border osize]);
      set(obj.gui.type_SEXT,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      obj.guiCreateCheckbox('type_MULT','MULT',1,'oppan3',[border border*3+osize*2 1-2*border osize]);
      set(obj.gui.type_MULT,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      obj.guiCreateCheckbox('type_COR','COR',1,'oppan3',[border border*4+osize*3 1-2*border osize]);
      set(obj.gui.type_COR,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      % - make selection
      osize=(1-3*border)/2;
      obj.guiCreatePushbutton('select_accept','Accept','oppan4',[border border 1-2*border osize]);
      set(obj.gui.select_accept,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      set(obj.gui.select_accept,'BackgroundColor','green');
      obj.guiCreatePushbutton('select_cancel','Cancel','oppan4',[border border*2+osize 1-2*border osize]);
      set(obj.gui.select_cancel,'Callback',@(src,event)guiIndexCallback(obj,src,event));
      set(obj.gui.select_cancel,'BackgroundColor','red');
      % fill gui
      obj.indexChoiceFromGui=obj.useCntrl;
      obj.indexChanChoiceFromGui=obj.useCntrlChan;
      obj.updateIndexGui;
    end
    function guiChanSelCallback(obj,src,~)
      sel=get(src,'Value');
      str=get(src,'String'); tstr=str{sel};
      csel=get(obj.gui.indexCbox1,'Value');
      mI=get(obj.gui.indexCbox1,'UserData');
      if ~get(obj.gui.showchan,'Value')
        tgl=true;
      else
        tgl=false;
      end
      for isel=csel
        if ismember(tstr,{'main' 'ampl'})
          obj.indexChanChoiceFromGui{mI(isel)}(1)=tgl;
        elseif strcmp(tstr,{'pha'})
          obj.indexChanChoiceFromGui{mI(isel)}(2)=tgl;
        else
          obj.indexChanChoiceFromGui{mI(isel)}(ismember({'x' 'dx' 'y' 'dy' 'z' 'dz'},tstr))=tgl;
        end
      end
      obj.updateIndexGui;
    end
    function guiIndexChanSelCallback(obj,src,~) %#ok<MANU>
      if get(src,'Value')
        set(src,'String','UnSelect Channel(s)')
      else
        set(src,'String','Select Channel(s)')
      end
    end
    function guiCselCallback(obj,~,~)
      sel=get(obj.gui.indexCbox1,'Value');
      mI=get(obj.gui.indexCbox1,'UserData');
      chnames={{'main'} {'x' 'dx' 'y' 'dy' 'z' 'dz'} {'ampl' 'pha'}};
      pslist=obj.PS_list; girlist=obj.GIRDER_list;
      allchans={};
      hwchans=obj.INDXchannels;
      for isel=sel
        for ichan=find(obj.useCntrlChan{mI(isel)})
          if ismember(ichan,hwchans{1,mI(isel)}) || ismember(ichan,hwchans{2,mI(isel)})
            if ismember(mI(isel),pslist)
              allchans{end+1}=chnames{1}{ichan};
            elseif ismember(mI(isel),girlist)
              allchans{end+1}=chnames{2}{ichan};
            else
              allchans{end+1}=chnames{3}{ichan};
            end
          end
        end
      end
      if isempty(allchans)
        set(obj.gui.chansel,'String','---')
      else
        set(obj.gui.chansel,'String',unique(allchans))
      end
    end
  end
end
