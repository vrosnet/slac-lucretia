classdef Track
% Lucretia Track Class: Main Lucretia tracking interface
  properties
    startInd
    finishInd
    BeamIn
    BeamOut
    loopFlag=0;
    firstBunch=1;
    lastBunch=1;
    instdata
    trackStatus
    statusCheck=true;
  end
  methods
    function obj=Track
      global BEAMLINE
      obj.startInd=1;
      obj.finishInd=length(BEAMLINE);
    end
    function obj=run(obj,BEAMLINE)
%       global BEAMLINE PS GIRDER KLYSTRON WF %#ok<*NUSED>
      % tracking in parallel across multiple hosts
      if strcmp(class(BEAMLINE),'Composite')
        startInd=obj.startInd;
        finishInd=obj.finishInd;
        BeamIn=obj.BeamIn;
        firstBunch=obj.firstBunch;
        lastBunch=obj.lastBunch;
        loopFlag=obj.loopFlag;
        statusCheck=obj.statusCheck;
        bb=BEAMLINE{1};
        spmd
          disp(isequal(bb,BEAMLINE{1}))
          [stat beamout]=TrackThru(1,100,BeamIn,1,1,0);
%           if statusCheck && trackStatus{1}~=1
%             error(trackStatus{2})
%           end
        end
        obj.BeamOut=BeamOut;
        obj.instdata=instdata;
        obj.trackStatus=trackStatus;
      else % local tracking
        [trackStatus BeamOut instdata]=TrackThru(obj.startInd,obj.finishInd,obj.BeamIn,obj.firstBunch,obj.lastBunch,obj.loopFlag);
        obj.trackStatus=trackStatus;
        obj.BeamOut=BeamOut;
        obj.instdata=instdata;
        if obj.statusCheck && trackStatus{1}~=1
          error(trackStatus{2})
        end
      end
    end
  end
  methods (Static)
    function trackVer=version
      trackVer=TrackThru('version');
    end
    function clear
      TrackThru('clear');
    end
  end
end