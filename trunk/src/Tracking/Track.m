classdef Track < handle
  % TRACK Lucretia beam tracking interface
  %   Perform tracking, either singly or on distributed Lucretia interface
  %   defined by optionally passed distributedLucretia object on construction
  %
  % Supports both asynchronous and synchronous parallel tracking (set in
  % distributedLucretia class passed to this object upon creation)
  %   - Asynchronous is slower to run than synchronous due to extra setup
  %   time but the trackThru command returns immediately which is useful if
  %   you want to process other commands serially whilst the parallel
  %   tracking is being computed
  %
  % Contructor:
  %   T=Track(InputBeam,distributedLucretiaObject)
  %     Omit distributedLucretiaObject if using this object in a non-parallel
  %     environment
  %
  % Main public methods (see doc help for details):
  %   trackThru - main tracking method
  %   trackThru('singleRay') - track single particle (mean of particle
  %                            distributions, sum of charge)
  %
  % Example:
  %  % Create a distributedLucretia object (choose synchronicity with
  %    isasyn=true|false)
  %  DL=distributedLucretia(isasyn)
  %  % Create Track object with a Lucretia InputBeam
  %  T=Track(InputBeam,DL) % Now set track indices, T.startInd,T.finishInd etc as desired
  %  % Make any lattice changes
  %  DL.latticeSyncVals(1).PS(53)=0.85;
  %  DL.latticeSyncVals(2).PS(53)=0.85;
  %  DL.PSTrim(53);
  %  % Issue track command (sends tracking job to parallel worker nodes)
  %  T.trackThru;
  %  % Wait for results (if isasyn=true this is instantaneous and command is
  %                      not necessary)
  %  DL.asynWait
  %  % Get results (the main output arguments from Lucretia's TrackThru
  %  %  function), if there were any tracking errors trying to access these
  %  %  parameters results in an error with the error output messages from the
  %  %  workers shown.
  %  for iw=DL.workers
  %    beamOut(iw)=T.beamOut{iw};
  %    trackStatus{iw}=T.trackStatus{iw};
  %    instrumentData(iw)=T.instrData{iw};
  %  end
  %
  % See also:
  %  TrackThru distributedLucretia
  %
  % Reference page in Help browser for list of accessible properties and
  % methods:
  %   <a href="matlab:doc Track">doc Track</a>
  %
  % Full lucretia documentation available online:
  %   <a href="http://www.slac.stanford.edu/accel/ilc/codes/Lucretia">Lucretia</a>
  properties
    startInd % Finish tracking index
    finishInd % Start tracking index
    firstBunch=1; % first bunch to track (if multibunch)
    lastBunch=1; % last bunch to track (if multibunch)
    loopFlag=0; % loop over elements (0) or over bunches (1)
    beamType=0; % 0=all input beams the same, 1=possibly different beams for each worker
    csrStoreData=false; % true: store Wakefield etc data at each CSR calculation point
    csrNbins=600; % number of histogram bins to use for CSR calculations
    verbose=0; % verbosity level (0= don't print anything, 1=print at each CSR integration step)
    csrData % Contains CSR data at each track point requiring CSR calculation if csrStoreData set true
    csrSmoothVal=3; % Amount of smoothing to apply to charge distribution for CSR calculation (integer >=1)
  end
  properties(SetAccess=protected)
    isDistrib % Is this Track object opererating in distributed mode?
    DL % distributedLucretia object reference
  end
  properties(Access=protected)
    dBeamIn % distributed input beam
    lBeamIn % local input beam
    instr
    beamout
    stat
    nray
    beamInSingle
  end
  properties(Dependent)
    beamIn % Lucretia beam structure to track
  end
  properties(Dependent,SetAccess=protected)
    instrData % (distributed) instrument data from tracking
    beamOut % (distributed) beam object post tracking
    trackStatus % (distributed) Lucretia status from tracking
  end
  
  %% Get/Set methods
  methods
    function set.trackStatus(obj,val)
      obj.stat=val;
    end
    function val=get.trackStatus(obj)
      if obj.isDistrib && ~obj.DL.synchronous
        val=asynGetData(obj.DL,1);
      else
        val=obj.stat;
      end
    end
    function set.instrData(obj,val)
      obj.instr=val;
    end
    function val=get.instrData(obj)
      if obj.isDistrib && ~obj.DL.synchronous
        val=asynGetData(obj.DL,3);
      else
        val=obj.instr;
      end
    end
    function set.beamOut(obj,val)
      obj.beamout=val;
    end
    function val=get.beamOut(obj)
      if obj.isDistrib && ~obj.DL.synchronous
        val=asynGetData(obj.DL,2);
      else
        val=obj.beamout;
      end
    end
    function set.beamIn(obj,beamStruc)
      % Update single ray structure
      if length(beamStruc.Bunch.Q)>1
        obj.beamInSingle=beamStruc;
        obj.beamInSingle.Bunch.Q=sum(beamStruc.Bunch.Q);
        obj.beamInSingle.Bunch.stop=0;
        obj.beamInSingle.Bunch.x=mean(beamStruc.Bunch.x,2);
      end
      % Deal with new Beam
      if obj.isDistrib && obj.DL.synchronous
        useWorkers=obj.DL.workers;
        spmd
          if ismember(labindex,useWorkers)
            BeamIn=beamStruc;
          end
        end
        obj.dBeamIn=BeamIn;
      elseif obj.isDistrib % asyn
        vars=whos('-file',obj.DL.asynDataFile);
        if any(ismember({vars.name},'trackBeamIn')) && obj.beamType
          ld=load(obj.DL.asynDataFile,'trackBeamIn');
          trackBeamIn=ld.trackBeamIn;
        end
        if obj.beamType
          for iw=obj.DL.workers
            trackBeamIn(iw)=beamStruc;
          end
        else
          trackBeamIn=beamStruc; %#ok<NASGU>
        end
        save(obj.DL.asynDataFile,'trackBeamIn','-append')
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
  
  %% Main public methods
  methods
    function obj=Track(beamIn,distribObj)
      global BEAMLINE
      if exist('distribObj','var') && ~isempty(distribObj)
        if ~strcmp(class(distribObj),'distributedLucretia') %#ok<STISA>
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
      obj.nray=numel(beamIn.Bunch.Q);
    end
    function trackThru(obj,cmd)
      % Asking for single-ray tracking?
      if exist('cmd','var') && isequal(cmd,'singleRay')
        doSingleRay=true;
      else
        doSingleRay=false;
      end
      if doSingleRay
        BeamIn=obj.beamInSingle;
      else
        BeamIn=obj.beamIn;
      end
      % tracking in parallel across possibly multiple workers
      if obj.isDistrib
        % If asynchronous, submit jobs and return
        if ~obj.DL.synchronous
          % Remove last job if there is one
          obj.DL.clearAsynJob;
          % Make new asyn job
          for iw=obj.DL.workers
            obj.DL.createAsynTask(@Track.asynTrack,3,{obj.DL.asynDataFile,iw,obj.startInd,obj.finishInd,...
              obj.firstBunch,obj.lastBunch,obj.loopFlag,doSingleRay});
          end
          obj.DL.launchAsynJob();
        else % if sycnchronous, submit tracking tasks to the pool and wait for them to finish
          startInd=obj.startInd;
          finishInd=obj.finishInd;
          b1=obj.firstBunch;
          b2=obj.lastBunch;
          lf=obj.loopFlag;
          verbose=obj.verbose;
          csrNbins=obj.csrNbins;
          csrSmoothVal=obj.csrSmoothVal;
          csrStoreData=obj.csrStoreData;
          useWorkers=obj.DL.workers;
          if ~isempty(indcsr) && obj.nray>1000
            spmd
              if ismember(labindex,useWorkers)
                [stat beamout instdata csrData]=Track.csrTrackThru(startInd,finishInd,BeamIn,b1,b2,...
                  lf,t0,indcsr,verbose,csrNbins,csrSmoothVal,csrStoreData,1) ;
              end
            end
            obj.csrData=csrData;
          else
            spmd
              if ismember(labindex,useWorkers)
                [stat beamout instdata]=TrackThru(startInd,finishInd,BeamIn,b1,b2,lf);
              end
            end
          end
          % Store results
          obj.trackStatus=stat;
          obj.beamOut=beamout;
          obj.instrData=instdata;
        end
      else % local tracking
        [stat beamout instdata]=TrackThru(obj.startInd,obj.finishInd,BeamIn,obj.firstBunch,obj.lastBunch,obj.loopFlag);
        obj.trackStatus=stat;
        obj.beamOut=beamout;
        obj.instrData=instdata;
      end
    end
  end
  
  %% Static methods (those needing to be called in worker environment)
  methods(Static)
    function [stat beamout instdata]=asynTrack(dataFile,iworker,i1,i2,b1,b2,lf,doSingleParticle)
      [BEAMLINE PS GIRDER KLYSTRON WF]=distributedLucretia.asynLoadLattice(dataFile,iworker); %#ok<NASGU,ASGLU>
      load(dataFile,'trackBeamIn')
      if length(trackBeamIn)>1
        beam=trackBeamIn(iworker);
      else
        beam=trackBeamIn;
      end
      if doSingleParticle
        beamInSingle=beam;
        beamInSingle.Bunch.Q=sum(beam.Bunch.Q);
        beamInSingle.Bunch.stop=0;
        beamInSingle.x=mean(beam.Bunch.x,2);
        beam=beamInSingle;
      end
      [stat beamout instdata]=TrackThru(i1,i2,beam,b1,b2,lf);
    end
  end
  
  methods(Access=private)
    function [bininds z Z ZSP]=doBinning(beamZ,nbin)
      zmin=min(-beamZ);
      zmax=max(-beamZ);
      if zmin==zmax
        error('Need some spread in z-distribution of bunch to compute CSR!')
      end
      z=linspace(zmin,zmax,nbin);
      z=z-mean(z);
      [~,bininds] = histc(-beamZ,z);
      [Z ZSP]=meshgrid(z,z);
    end
  end
  
end