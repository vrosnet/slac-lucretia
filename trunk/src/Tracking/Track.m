classdef Track < handle
% TRACK Lucretia beam tracking interface
%   Perform tracking, either singly or on distributed Lucretia interface
%   defined by optionally passed distributedLucretia object on construction
  properties
    startInd % Finish tracking index
    finishInd % Start tracking index
    firstBunch=1; % first bunch to track (if multibunch)
    lastBunch=1; % last buncg to track (if multibuncgh
    loopFlag=0; % loop over elements (0) or over bunches (1)
    beamType=0; % 0=all input beams the same, 1=possibly different beams for each worker
  end
  properties(SetAccess=private)
    instrData % (distributed) instrument data from tracking
    trackStatus % (distributed) Lucretia status from tracking
    beamOut % (distributed) beam object post tracking
    isDistrib % Is this Track object opererating in distributed mode?
    DL % distributedLucretia object reference
  end
  properties(Access=private)
    dBeamIn % distributed input beam
    lBeamIn % local input beam
    trackJob % distributed job data
  end
  properties(Dependent)
    beamIn % Lucretia beam structure to track
  end
  properties(Dependent,SetAccess=private)
    asynJobStatus='pending' % Job status for asyn jobs= 'pending','runnning','queued','finished','failed','destroyed','unavailable'
  end
  
  % Get/Set methods
  methods
    function status=get.asynJobStatus(obj)
      try
        status=obj.trackJob.State;
      catch
        status='unknown';
      end
    end
    function set.beamIn(obj,beamStruc)
      if obj.isDistrib && obj.DL.synchronous
        useWorkers=obj.DL.workers;
        spmd
          if ismember(labindex,useWorkers)
            BeamIn=beamStruc;
          end
        end
        obj.dBeamIn=BeamIn;
      elseif obj.isDistrib % asyn
        dataFile=fullfile(obj.DL.sched.DataLocation,'distributedLucretiaStartupData.mat');
        vars=whos('-file',dataFile);
        if any(ismember({vars.name},'trackBeamIn')) && obj.beamType
          ld=load(dataFile,'trackBeamIn');
          trackBeamIn=ld.trackBeamIn;
        end
        if obj.beamType
          for iw=obj.DL.workers
            trackBeamIn(iw)=beamStruc;
          end
        else
          trackBeamIn=beamStruc; %#ok<NASGU>
        end
        save(dataFile,'trackBeamIn','-append')
      else
        if obj.beamType
          obj.lBeamIn=[];
          obj.lBeamIn(obj.DL.workers)=beamStruc;
        else
          obj.lBeamIn=beamStruc;
        end
      end
    end
    function beam=get.beamIn(obj)
      if obj.isDistrib && obj.DL.synchronous
        beam=obj.dBeamIn;
      elseif obj.isDistrib
        dataFile=fullfile(obj.DL.sched.DataLocation,'distributedLucretiaStartupData.mat');
        vars=whos('-file',dataFile);
        if any(ismember({vars.name},'trackBeamIn'))
          ld=load(dataFile,'trackBeamIn');
          beam=ld.trackBeamIn;
        else
          beam=[];
        end
      else
        beam=obj.lBeamIn;
      end
    end
  end
  
  % Main public methods
  methods
    function obj=Track(beamIn,distribObj)
      global BEAMLINE
      if exist('distribObj','var')
        if ~strcmp(class(distribObj),'distributedLucretia')
          error('Can only pass a distributedLucretia object to Track')
        end
        obj.isDistrib=true;
        obj.DL=distribObj;
      else
        obj.isDistrib=false;
      end
      if ~exist('beamIn','var')
        error('Must supply input beam structure')
      end
      obj.startInd=1;
      obj.finishInd=length(BEAMLINE);
      obj.beamIn=beamIn;
    end
    function asynWait(obj,timeout)
      % asynWait(obj,timeout)
      %  Wait for asynchronous job to complete
      %  Optional timeout in seconds
      if obj.DL.synchronous
        return
      end
      t0=clock;
      while ~strcmp(obj.asynJobStatus,'finished') && ~strcmp(obj.asynJobStatus,'failed') && ~strcmp(obj.asynJobStatus,'destroyed') ...
          && ~strcmp(obj.asynJobStatus,'unavailable')
        pause(1)
        if exist('timeout','var') && etime(clock,t0)>timeout
          break
        else
          t0=clock;
        end
      end
    end
    function trackThru(obj)
      % tracking in parallel across possibly multiple workers
      if obj.isDistrib
        % Remove last job if there is one
        if ~isempty(obj.trackJob)
          try
            obj.trackJob.destroy;
          catch
          end
        end
        % If asynchronous, submit jobs and return
        if ~obj.DL.synchronous
          obj.trackJob=createJob(obj.DL.sched);
          dataFile=fullfile(obj.DL.sched.DataLocation,'distributedLucretiaStartupData.mat');
          for iw=obj.DL.workers
            obj.trackJob.createTask(@asynTrack,4,{dataFile,iw,obj.startInd,obj.finishInd,...
              obj.firstBunch,obj.lastBunch,obj.loopFlag,obj.DL.latticeSyncVals(iw)});
          end
          obj.trackJob.submit;
        else % if sycnchronous, submit tracking tasks to the pool and wait for them to finish
          startInd=obj.startInd;
          finishInd=obj.finishInd;
          BeamIn=obj.beamIn;
          b1=obj.firstBunch;
          b2=obj.lastBunch;
          lf=obj.loopFlag;
          useWorkers=obj.DL.workers;
          spmd
            if ismember(labindex,useWorkers)
              [stat beamout instdata]=TrackThru(startInd,finishInd,BeamIn,b1,b2,lf);
            end
          end
          % Store results
          obj.trackStatus=stat;
          obj.beamOut=beamout;
          obj.instrData=instdata;
        end 
      else % local tracking
        [stat beamout instdata]=TrackThru(obj.startInd,obj.finishInd,obj.beamIn,obj.firstBunch,obj.lastBunch,obj.loopFlag);
        obj.trackStatus=stat;
        obj.beamOut=beamout;
        obj.instrData=instdata;
      end
    end
  end
  
  % Static methods (those needing to be called in worker environment)
  methods(Static)
    function [stat beamout instdata retVals]=asynTrack(dataFile,iworker,i1,i2,b1,b2,lf,syncVals,method)
      load(dataFile,'BEAMLINE','PS','GIRDER','WF','KLYSTRON','trackBeamIn');
      [PS,KLYSTRON,GIRDER,retVals]=distributedLucretia.syncLattice(syncVals,PS,KLYSTRON,GIRDER,method); %#ok<NODEF>
      if length(trackBeamIn)>1
        beam=trackBeamIn(iworker);
      else
        beam=trackBeamIn;
      end
      [stat beamout instdata]=TrackThru(i1,i2,beam,b1,b2,lf);
    end
  end
end