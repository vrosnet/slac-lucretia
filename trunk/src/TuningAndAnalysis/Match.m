classdef Match
  %MATCH Lucretia beam matching using GetTwiss and TrackThru
  %
  %Useage example below (see 'doc Match' for complete list of properties and methods):
  %Create object:
  %M=Match;
  %
  %Assign Lucretia beam if tracking needed (matching S,T or U):
  %M.beam=Beam;
  %
  %Assign Lucretia Initial structure if twiss matching needed
  % (anything other than S,T or U):
  %  M.initStruc=InitialStructure;
  %  
  %Assign initial track/twiss point corresponding to above:
  %M.iInitial=findcells(BEAMLINE,'Name','InitialMarkerName');
  %
  %Choose an optimizer to use:
  %M.optim='lsqnonlin';
  %  Supported optimizers are:
  %    * fminsearch (supported in standard Matlab)
  %    * fmincon (requires optimization toolbox)
  %    * lsqnonlin (requries optimization toolbox)
  %    * fgoalattain (requires optimization toolbox)
  %    * genetic (requires global optimization and optimization toolboxes)
  %  The recommended optimizers are lsqnonlin for twiss optimization only
  %  and fminsearch or fmincon for tracking-based optimization (S, T or U).
  %  The Genetic algorithm requires a lot of tuning of its input parameters
  %  to give sensible answers, you will need to play with the gaoptimset options
  %  inside the doMatch method.
  %  
  %Add optimization variables
  %M=addVariable(M,'PS',BEAMLINE{quadEle}.PS,'Ampl',0.8,1.2);
  %  Available variable types are: PS,BEAMLINE,GIRDER,KLYSTRON
  %  You can use any available field of the supported types that contains a
  %  scalar quantity. Repeat above for each required variable. The last 2
  %  input parameters define the lower and upper bounds for this variable.
  %  For the constrained optimizers the final variable values are
  %  guaranteed to be within this range, for the unconstrained ones this
  %  range is used as part of the weighting system to restrain the allowed
  %  movement of this variable.
  %
  %Add constraints
  %M=addMatch(M,beamlineElement1,'alpha_y',0,1e-4);
  %M=addMatch(M,beamlineElement2,'Sigma',35e-9^2,1e-9^2,'33');
  %  Add a constraint for the optimizer at the given BEAMLINE element
  %  number. Available constraint types are: alpha_x, beta_x, eta_x,
  %  etap_x, nu_x, NEmit_x, S, T, U. Also _y for vertical twiss parameters.
  %  For S, T, U (second moment, and second and third-order beam
  %  correlations) also supply correlation matrix elements as last
  %  argument. e.g. '33' for sigma_33. '3126' for U3126 etc... Third
  %  argument is required match value. Forth argument is tolerance, the
  %  matcher will stop when all added constraints are closer than this
  %  value to the desired match value.
  %  Repeat for all required match parameters
  %
  %display(M) [or simply type 'M']
  %  Look at the list of assigned variables and constraints and the current
  %  match data
  %
  %M=M.doMatch;
  %  Perform the matching
  %
  %display(M)
  %  See how it did
  %
  % (Note: to use the T/U match conditions you need the "polyfitn" tools
  % supplied in the Lucretia/src/utils/PolyfitnTools directory)
  %
  % See also:
  %  InitCondStruc MakeBeam6dGauss findcells fminsearch fmincon lsqnonlin
  %  fgoalattain
  %
  % Reference page in Help browser for list of accessible properties and
  % methods:
  %   <a href="matlab:doc Match">doc Match</a>
  
  properties
    beam % Beam used for tracking (should be Lucretia Beam format used by TrackThru)
    initStruc % Lucretia Initial structure at iInitial BEAMLINE location
    iInitial % Initial beamline index (corresponding to beam)
    optim='lsqnonlin'; % Optimizer to use (options are fminsearch, fmincon, lsqnonlin, genetic, fgoalattain)
    verbose=false; % if set true, then display match conditions as optimiser runs
    optimDisplayRate=10; % How often to update the status display if verbose set to true (seconds)
    optimDisplay='final'; % what to pass to the optimiser display parameter (e.g. 'off', 'final', 'iter')
    useParallel=false; % Use parallel processing for optimisation. Requires parallel toolbox and fmincon or genetic or fgoalattain algrithm use.
    storeOptim=false; % Store all optimisation step data
  end
  properties(SetAccess=private)
    matchType={}; % Cell array of match type descriptors (see allowedMatchTypes property)
    matchTypeQualifiers={}; % Qualifiers for given types of match (e.g. '13' tp set Sigma_13, '166' to set T_166), cell array matching matchType size
    matchWeights=[]; % Vector of weights for match entries
    matchVals % The values for the match constraints specified
    matchInd % Beamline vector of indices defining match points
    varType={}; % Cell array of variable types (see allowedVariableTypes property)
    varIndex={}; % Cell array of indexes pertaining to varType (if cell entry is vector, then tie together multiple elements using same varField)
    varField={}; % Field of requested variable (must be existing field of requested variable type), cell array of same dimension as varType
    varFieldIndex % Index of field
    varLimits % array of high/low limits for variables [2,length(obj.varType)]
  end
  properties(Dependent)
    dotrack
    dotwiss
    varVals
    iMatch % unique list of required match points
    optimVals % match values returned from the optimizer
    Ix % horizontal Initial vals for GetTwiss
    Iy % vertical Initial vals for GetTwiss
    optimData % contraint, variable and optimStop values at each iteration step
  end
  properties(Constant)
    allowedMatchTypes={'alpha_x' 'alpha_y' 'beta_x' 'beta_y' 'NEmit_x' 'NEmit_y' 'eta_x' 'etap_x' 'eta_y' 'etap_y' 'nu_x' 'nu_y' 'Sigma' 'T' 'U'};
    allowedVariableTypes={'PS' 'GIRDER' 'KLYSTRON' 'BEAMLINE'};
  end
  
  methods
    function obj=Match
      % Match - checks for optimization toolbox upon object creation
      global MATCHDATA
      vout=ver;
      if any(ismember({vout(:).Name},'Optimization Toolbox'))
        obj.optim='lsqnonlin';
      else
        obj.optim='fminsearch';
      end
      MATCHDATA=[];
    end
    function obj=set.useParallel(obj,val)
      if val
        v=ver;
        if ~any(arrayfun(@(x) strcmp(v(x).Name,'Parallel Computing Toolbox'),1:length(v)))
          error('Must have parallel computing toolbox to use parallel option')
        end
        if ~ismember(obj.optim,{'fmincon','genetic','fgoalattain'}) %#ok<MCSUP>
          error('Must use one of the algorithms: ''fmincon'',''genetic'',''fgoalattain'' to use parallel option')
        end
      end
      obj.useParallel=val;
    end
    function obj=addMatch(obj,mind,type,vals,tol,typeQual)
      % ADDMATCH - add a match condition to the problem set
      % addMatch(obj,type,weight,vals,tol [,typeQual])
      % type = string or cell list of strings, see allowedMatchTypes property for available list
      % weight = double, weighting function for constraints
      % vals = double, desired match values
      % typeQual = qualifier (optional) string or cell vector of strings
      %            (e.g. '13' tp set Sigma_13, '166' to set T_166)
      % tol = tolerance for this match condition
      % ... All above can be scalar or vectors of same length to either add
      %     a single match condition or multiple
      if ~exist('type','var') || any(~ismember(type,obj.allowedMatchTypes))
        error('Must choose from available list of match types')
      end
      if ~exist('mind','var') || (iscell(type) && length(mind)~=length(type))
        error('Must supply BEAMLINE index or vector of indices for this match constraint (length same as number of type declarations)')
      end
      if ~exist('vals','var') || ~isnumeric(vals) || (iscell(type) && length(vals)~=length(type))
        error('Must supply match value(s) ''var'' scalar double or vector of same length as type')
      end
      if ~exist('tol','var') || ~isnumeric(tol) || (iscell(type) && length(tol)~=length(type))
        error('Must supply match value(s) ''tol'' scalar double or vector of same length as type')
      end
      if iscell(type)
        nadd=length(type);
      else
        nadd=1;
      end
      if exist('typeQual','var') && iscell(typeQual) && iscell(type) && length(typeQual) ~= length(type)
        error('If supply typeQual, must be either scalar cell or vector same length as ''type''')
      elseif exist('typeQual','var') && iscell(typeQual)
        obj.matchTypeQualifiers(end+1:end+nadd)=typeQual;
      elseif exist('typeQual','var')
        obj.matchTypeQualifiers{end+1}=typeQual;
      elseif ~exist('typeQual','var')
        for iadd=1:nadd
          obj.matchTypeQualifiers{end+1}='';
        end
      end
      if iscell(type)
        obj.matchType(end+1:end+nadd)=type;
      else
        obj.matchType{end+1}=type;
      end
      obj.matchWeights(end+1:end+nadd)=tol;
      obj.matchVals(end+1:end+nadd)=vals;
      obj.matchInd(end+1:end+nadd)=mind;
    end
    function obj=addVariable(obj,type,typeIndex,field,lowLimit,highLimit,fieldIndex)
      % ADDVARIABLE - add a variable to match problem set
      % addVariable(obj,type,typeIndex,field,lowLimit,highLimit [,fieldIndex])
      % type = string, or cell vector of strings: see allowedVariableTypes property for allowed options
      % typeIndex = integer
      % field = string, or cell of strings length of type: fieldname of type
      % lowLimit = lower bound on possible vals
      % highLimit = upper bound on possible vals
      % fieldIndex (optional) = double, index of field
      % ... All above can be scalar or vectors of same length for single or
      %     multiple variable declarations
      if ~exist('type','var') || any(~ismember(type,obj.allowedVariableTypes))
        error('Must supply variable type(s) from allowed list')
      end
      if ~exist('typeIndex','var') || (iscell(type) && length(typeIndex)~=length(type))
        error('Must supply typeIndex integer or vector of integers same length as type')
      end
      if ~exist('field','var')  || (iscell(type) && length(field)~=length(type))
        error('Must supply field name or cell of field names same length as type')
      end
      if ~exist('lowLimit','var') || (iscell(type) && length(lowLimit)~=length(type))
        error('Must supply low limit to variable or vector of the same length as type')
      end
      if ~exist('highLimit','var') || (iscell(type) && length(highLimit)~=length(type))
        error('Must supply high limit to variable or vector of the same length as type')
      end
      if iscell(type)
        nadd=length(type);
      else
        nadd=1;
      end
      if exist('fieldIndex','var') && iscell(type) && length(fieldIndex)~=length(type)
        error('If supplying fieldIndex, must be of same length as number of types supplied')
      elseif exist('fieldIndex','var')
        obj.varFieldIndex(end+1:end+nadd)=fieldIndex;
      elseif ~exist('fieldIndex','var')
        obj.varFieldIndex(end+1:end+nadd)=1;
      end
      if iscell(type)
        obj.varType(end+1:end+nadd)=type;
      else
        obj.varType{end+1}=type;
      end
      if iscell(typeIndex)
        obj.varIndex(end+1:end+nadd)=typeIndex;
      else
        obj.varIndex{end+1}=typeIndex;
      end
      if iscell(field)
        obj.varField(end+1:end+nadd)=field;
      else
        obj.varField{end+1}=field;
      end
      obj.varLimits(1,end+1:end+nadd)=lowLimit;
      obj.varLimits(2,end-nadd+1:end)=highLimit;
    end
    function display(obj)
      fprintf('Lucretia match using optimizer: %s\n',obj.optim)
      fprintf('-----------------------------------------\n')
      fprintf('Beam sizes at match points:\n')
      ov=obj.optimVals;
      v=obj.minFunc('getvals');
      for im=1:length(obj.iMatch)
        fprintf('%d: sigma_x= %g sigma_y = %g\n',obj.iMatch(im),v.bs.x(im),v.bs.y(im))
      end
      fprintf('-----------------------------------------\n')
      fprintf('Constraints / Desired Vals / Current Vals / Weights\n')
      fprintf('\n')
      for itype=1:length(obj.matchType)
        fprintf('%d: %s (%s) / %g / %g / %g\n',obj.matchInd(itype),obj.matchType{itype},obj.matchTypeQualifiers{itype},...
          obj.matchVals(itype),ov(itype),obj.matchWeights(itype))
      end
      fprintf('\n')
      fprintf('-----------------------------------------\n')
      fprintf('Variables / Val / Low / High\n')
      fprintf('\n')
      curVals=obj.varVals;
      for itype=1:length(obj.varType)
        indstr=' ';
        for ind=obj.varIndex{itype}
          indstr=[indstr num2str(ind) ' '];
        end
        fprintf('%s(%s).%s(%d) / %g / %g / %g\n',obj.varType{itype},indstr,obj.varField{itype},...
          obj.varFieldIndex(itype),curVals(itype),obj.varLimits(1,itype),obj.varLimits(2,itype))
      end
      fprintf('============================================\n')
      fprintf('============================================\n')
    end
    function obj=doMatch(obj)
      % Perform programmed match procedure
      global BEAMLINE PS GIRDER KLYSTRON
      % Check structure good for fitting
      if obj.dotrack && isempty(obj.beam)
        error('Must supply Lucretia beam for tracking')
      end
      if obj.dotwiss && isempty(obj.initStruc)
        error('Must supply Lucretia Initial structure for Twiss propagation')
      end
      % --- Perform matching
      % Set global options
      if strcmp(obj.optim,'genetic')
        opts=gaoptimset('Display',obj.optimDisplay,...
          'TolCon',1e-6,'TolFun',1e-6,'Generations',100000,'PopulationSize',100,'UseParallel','never','Vectorized','off');
      elseif strcmp(obj.optim,'lsqnonlin')
        opts=optimset('Display',obj.optimDisplay,'MaxFunEvals',200000,'MaxIter',100000,'TolX',1e-6,'TolFun',1e-6,...
          'OutputFcn',@(x,optimValues,state) optimOutFun(obj,x,optimValues,state));
      elseif strcmp(obj.optim,'fmincon')
        opts=optimset('Display',obj.optimDisplay,'OutputFcn',@(x,optimValues,state) optimOutFun(obj,x,optimValues,state),'MaxFunEvals',100000,...
          'MaxIter',100000,'TolX',1e-6,'TolFun',1,'UseParallel','never','Algorithm','active-set');
      else
        opts=optimset('Display',obj.optimDisplay,'OutputFcn',@(x,optimValues,state) optimOutFun(obj,x,optimValues,state),'MaxFunEvals',100000,...
          'MaxIter',100000,'TolX',1e-6,'TolFun',1,'UseParallel','never');
      end
      % Setup parallel envirnoment
      varVals=obj.varVals;
      if obj.useParallel
        save blvals BEAMLINE PS GIRDER KLYSTRON
        matlabpool open
        spmd
          bl=load('blvals');
          BEAMLINE=bl.BEAMLINE; %#ok<*SPGP>
          PS=bl.PS;
          GIRDER=bl.GIRDER;
          KLYSTRON=bl.KLYSTRON;
        end
      end
      
      % Perform fit
      % - normalise variables
      varW=(obj.varLimits(2,:)-obj.varLimits(1,:));
      switch obj.optim
        case 'genetic'
          x=ga(@(x) minFunc(obj,x,obj.dotrack,obj.dotwiss),length(varVals),[],[],[],[],...
            obj.varLimits(1,:)./varW,obj.varLimits(2,:)./varW,[],opts);
        case 'lsqnonlin'
          x=lsqnonlin(@(x) minFunc(obj,x,obj.dotrack,obj.dotwiss,varW,varVals),...
            varVals./varW,obj.varLimits(1,:)./varW,obj.varLimits(2,:)./varW,opts);
        case 'fminsearch'
          x=fminsearch(@(x) minFunc(obj,x,obj.dotrack,obj.dotwiss),varVals./varW,opts);
        case 'fmincon'
          x=fmincon(@(x) minFunc(obj,x,obj.dotrack,obj.dotwiss),varVals./varW,[],[],[],[],...
            obj.varLimits(1,:)./varW,obj.varLimits(2,:)./varW,[],opts);
        case 'fgoalattain'
          x=fgoalattain(@(x) minFunc(obj,x,obj.dotrack,obj.dotwiss),varVals,obj.matchVals,obj.matchWeights,...
            [],[],[],[],obj.varLimits(1,:),obj.varLimits(2,:),[],opts);
        otherwise
          error('Unknown or unsupported optimizer');
      end
      if obj.useParallel
        matlabpool close
        clear BEAMLINE PS GIRDER KLYSTRON
      end
      if strcmp(obj.optim,'fgoalattain')
        obj.varVals=x;
      else
        obj.varVals=x.*varW;
      end
    end
    function out=get.iMatch(obj)
      out=unique(obj.matchInd);
    end
    function out=get.optimData(obj)
      out=optimOutFun(obj,'getdata');
    end
    function out=get.dotrack(obj)
      % Check for match requirements to see if tracking is required
      if any(ismember(obj.matchType,{'Sigma' 'T' 'U' 'NEmit_x' 'NEmit_y'}))
        out=true;
      else
        out=false;
      end
    end
    function out=get.dotwiss(obj)
      % Check for match requirements to see if tracking is required
      if any(ismember(obj.matchType,{'alpha_x' 'alpha_y' 'beta_x' 'beta_y' 'eta_x' 'etap_x' 'eta_y' 'etap_y' 'nu_x' 'nu_y'}))
        out=true;
      else
        out=false;
      end
    end
    function vals=get.optimVals(obj)
      varW=(obj.varLimits(2,:)-obj.varLimits(1,:));
      if strcmp(obj.optim,'fgoalattain')
        minFunc(obj,obj.varVals,obj.dotrack,obj.dotwiss);
      else
        minFunc(obj,obj.varVals./varW,obj.dotrack,obj.dotwiss);
      end
      v=obj.minFunc('getvals');
      vals=v.func;
    end
    function vals=get.varVals(obj)
      global BEAMLINE PS GIRDER KLYSTRON
      for itype=1:length(obj.varType)
        ind=obj.varIndex{itype}(1);
        switch obj.varType{itype}
          case 'BEAMLINE'
            vals(itype)=BEAMLINE{ind}.(obj.varField{itype});
          case 'PS'
            vals(itype)=PS(ind).(obj.varField{itype});
          case 'GIRDER'
            vals(itype)=GIRDER{ind}.(obj.varField{itype});
          case 'KLYSTRON'
            vals(itype)=KLYSTRON(ind).(obj.varField{itype});
        end
      end
    end
    function obj=set.varVals(obj,x)
      % Set variables
      global BEAMLINE PS GIRDER KLYSTRON
      for itype=1:length(obj.varType)
        for ind=obj.varIndex{itype}
          switch obj.varType{itype}
            case 'BEAMLINE'
              BEAMLINE{ind}.(obj.varField{itype})=x(itype);
            case 'PS'
              PS(ind).(obj.varField{itype})=x(itype);
            case 'GIRDER'
              GIRDER{ind}.(obj.varField{itype})=x(itype);
            case 'KLYSTRON'
              KLYSTRON(ind).(obj.varField{itype})=x(itype);
          end
        end
      end
    end
  end
  methods(Access=private)
    function [stop,options,optchanged]=gaOptimOutFun(obj,options,~,~,~)
      optchanged=false;
      stop = optimOutFun(obj);
    end
    function stop = optimOutFun(obj,x,optimval,~)
      persistent lastUpdate
      global MATCHDATA
      stop=false;
      if (isfield(optimval,'fval') && optimval.fval<1) || (isfield(optimval,'resnorm') && optimval.resnorm<1)
        stop=true;
      end
      if obj.verbose && (isempty(lastUpdate) || etime(clock,lastUpdate)>obj.optimDisplayRate)
        display(obj);
        lastUpdate=clock;
      end
      if obj.storeOptim
        MATCHDATA(end+1).optimvals=vals.func;
        MATCHDATA(end).varvals=x;
      end
    end
    function F = minFunc(obj,x,dotrack,dotwiss,varW,varVals)
      % The minimizer function
      persistent lastvals bs
      
      % just get last match values
      if isequal(x,'getvals')
        F.func=lastvals;
        F.bs=bs;
        return
      end
      
      % Undo variable normalisation
      if ~strcmp(obj.optim,'fgoalattain')
        vm=(obj.varLimits(2,:)-obj.varLimits(1,:));
        x=x.*vm;
      end
      
      % Set the variables
      if obj.useParallel
        for itype=1:length(obj.varType)
          for ind=obj.varIndex{itype}
            switch obj.varType{itype}
              case 'BEAMLINE'
                BEAMLINE{ind}.(obj.varField{itype})=x(itype); %#ok<NASGU>
              case 'PS'
                PS(ind).(obj.varField{itype})=x(itype); %#ok<NASGU>
              case 'GIRDER'
                GIRDER{ind}.(obj.varField{itype})=x(itype); %#ok<NASGU>
              case 'KLYSTRON'
                KLYSTRON(ind).(obj.varField{itype})=x(itype); %#ok<NASGU>
            end
          end
        end
      else
        obj.varVals=x;
      end
      
      % Get twiss in correct format
      if dotwiss
        Ix.beta=obj.initStruc.x.Twiss.beta;
        Ix.alpha=obj.initStruc.x.Twiss.alpha;
        Ix.eta=obj.initStruc.x.Twiss.eta;
        Ix.etap=obj.initStruc.x.Twiss.etap;
        Ix.nu=obj.initStruc.x.Twiss.nu;
        Iy.beta=obj.initStruc.y.Twiss.beta;
        Iy.alpha=obj.initStruc.y.Twiss.alpha;
        Iy.eta=obj.initStruc.y.Twiss.eta;
        Iy.etap=obj.initStruc.y.Twiss.etap;
        Iy.nu=obj.initStruc.y.Twiss.nu;
      end
      
      % Get the data (do the tracking)
      F=zeros(1,length(obj.matchType));
      for itrack=1:length(obj.iMatch)
        try
          if itrack==1
            if dotrack
              [stat beamout]=TrackThru(obj.iInitial,obj.iMatch(1),obj.beam,1,1,0);
              if stat{1}~=1; error('Error in tracking in minimiser function: %s',stat{2}); end;
            end
            if dotwiss
              [stat T]=GetTwiss(obj.iInitial,obj.iMatch(1),Ix,Iy);
              if stat{1}~=1; error('Error in tracking in minimiser function: %s',stat{2}); end;
            end
          else
            if dotrack
              [stat beamout]=TrackThru(obj.iMatch(itrack-1),obj.iMatch(itrack),beamout,1,1,0);
              if stat{1}~=1; error('Error in tracking in minimiser function: %s',stat{2}); end;
            end
            if dotwiss
              %             Ix.beta=T.betax; Ix.alpha=T.alphax; Ix.eta=T.etax; Ix.etap=T.etapx; Ix.nu=T.nux;
              %             Iy.beta=T.betay; Iy.alpha=T.alphay; Iy.eta=T.etay; Iy.etap=T.etapy; Iy.nu=T.nuy;
              [stat T]=GetTwiss(obj.iInitial,obj.iMatch(itrack),Ix,Iy);
              if stat{1}~=1; error('Error in tracking in minimiser function: %s',stat{2}); end;
            end
          end
        catch
          F=ones(size(F)).*1e10;
          continue
        end
        % Extract the data
        for itype=1:length(obj.matchType)
          if obj.matchInd(itype)~=obj.iMatch(itrack); continue; end;
          switch obj.matchType{itype}
            case 'alpha_x'
              F(itype)=T.alphax(end);
            case 'alpha_y'
              F(itype)=T.alphay(end);
            case 'beta_x'
              F(itype)=T.betax(end);
            case 'beta_y'
              F(itype)=T.betay(end);
            case 'eta_x'
              F(itype)=T.etax(end);
            case 'eta_y'
              F(itype)=T.etay(end);
            case 'etap_x'
              F(itype)=T.etapx(end);
            case 'etap_y'
              F(itype)=T.etapy(end);
            case 'nu_x'
              F(itype)=T.nux(end);
            case 'nu_y'
              F(itype)=T.nuy(end);
            case 'NEmit_x'
              if ~exist('nx','var')
                [nx,ny] = GetNEmitFromBeam( beamout ,1);
              end
              F(itype)=nx;
            case 'NEmit_y'
              if ~exist('ny','var')
                [nx,ny] = GetNEmitFromBeam( beamout ,1);
              end
              F(itype)=ny;
            case 'Sigma'
              if ~exist('S','var')
                S=cov(beamout.Bunch.x');
              end
              F(itype)=S(str2double(obj.matchTypeQualifiers{itype}(1)),str2double(obj.matchTypeQualifiers{itype}(2)));
            case {'T' 'U'}
              if ~exist('fitCoef','var') || dim~=str2double(obj.matchTypeQualifiers{itype}(1))
                dim=str2double(obj.matchTypeQualifiers{itype}(1));
                [fitTerm,fitCoef] = beamTerms(dim,beamout);
              end
              term=[0 0 0 0]; termind=[1 2 3 4 6]; termind=termind(termind~=dim);
              for iterm=2:length(obj.matchTypeQualifiers{itype})
                term(termind==str2double(obj.matchTypeQualifiers{itype}(iterm)))=...
                  term(termind==str2double(obj.matchTypeQualifiers{itype}(iterm)))+1;
              end
              F(itype)=fitCoef(arrayfun(@(x) isequal(fitTerm(x,:),term),1:length(fitTerm)));
          end
        end
        
        % Beamsize data
        if dotrack
          bs.x(itrack)=std(beamout.Bunch.x(1,:));
          bs.y(itrack)=std(beamout.Bunch.x(3,:));
        else
          bs.x(itrack)=sqrt((obj.initStruc.x.NEmit/(obj.initStruc.Momentum/0.511e-3))*T.betax(end));
          bs.y(itrack)=sqrt((obj.initStruc.y.NEmit/(obj.initStruc.Momentum/0.511e-3))*T.betay(end));
        end
        
      end
      
      % Subtract desired values so minimiser does the correct thing
      % And apply weights
      lastvals=F;
      
      if strcmp(obj.optim,'lsqnonlin')
        F=(F-obj.matchVals)./obj.matchWeights;
        if exist('varW','var')
          F(end+1:end+length(varW))=(x-varVals)./varW;
        end
      elseif ~strcmp(obj.optim,'fgoalattain')
        F=sum(((F-obj.matchVals)./obj.matchWeights).^2);
      end 
    end
  end
  
end

