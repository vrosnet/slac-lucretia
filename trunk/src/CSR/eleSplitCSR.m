function eleSplitCSR(bsplitIn,dsplit,bunchLength)
%
global BEAMLINE


iele=1;
while iele<=length(BEAMLINE)
  
  lastEle=iele;
  if strcmp(BEAMLINE{iele}.Class,'SBEN')
    
    % Is this a split-element bend? (find end of physical magnet)
    magInd=iele;
    if iele<length(BEAMLINE)
      for ind=iele+1:length(BEAMLINE)
        if ind==length(BEAMLINE) || (~strcmp(BEAMLINE{ind}.Class,'SBEN') && isfield(BEAMLINE{ind},'L') && BEAMLINE{ind}.L>0)
          break
        elseif strcmp(BEAMLINE{ind}.Class,'SBEN')
          magInd(end+1)=ind;
        end
      end
    end
    
    % Split up all bends elements part of this physical magnet
    bsplit=floor(bsplitIn/length(magInd));
    magid=0; Angle=0; L=0;
    for imag=1:length(magInd)
      B=BEAMLINE{magInd(imag)+magid};
      Angle=Angle+B.Angle;
      L=L+B.L;
      Bs=B;
      Bs.L=B.L/bsplit;
      Bs.B=B.B./bsplit;
      Bs.Angle=B.Angle./bsplit;
      Bs.EdgeAngle=[0 0];
      Bs.FINT=[0 0];
      Bs.EdgeCurvature=[0 0];
      Bs.TrackFlag.doCSR=true;
      B1=Bs;
      B1.FINT=[B.FINT(1) 0];
      B1.EdgeAngle=[B.EdgeAngle(1) 0];
      B1.EdgeCurvature=[B.EdgeCurvature(1) 0];
      B1.TrackFlag.doCSR=true;
      B2=Bs;
      B2.FINT=[0 B.FINT(2)];
      B2.EdgeAngle=[0 B.EdgeAngle(2)];
      B2.EdgeCurvature=[0 B.EdgeCurvature(2)];
      B2.TrackFlag.doCSR=true;
      newBL={};
      newBL{end+1}=B1;
      for isplit=2:bsplit-1
        newBL{end+1}=Bs;
      end
      newBL{end+1}=B2;
      if magInd(imag)>1 && (magInd(imag)+1+magid)<length(BEAMLINE)
        BEAMLINE=[BEAMLINE(1:magInd(imag)+magid-1); newBL'; BEAMLINE(magInd(imag)+magid+1:end)];
      elseif magInd(imag)==1
        BEAMLINE=[newBL'; BEAMLINE(2:end)];
      elseif (magInd(imag)+1+magid)==length(BEAMLINE)
        BEAMLINE=[BEAMLINE(1:magInd(imag)+magid-1); newBL'];
      end
      magid=magid+length(newBL)-1;
    end
    SetSPositions( 1, length(BEAMLINE), 0 );
    lastEle=magInd(imag)+magid;
    
    % How far downstream to consider CSR effects?
    R=L/(2*sin(Angle/2));
    lDecay=3*(24*bunchLength*R^2)^(1/3);
    
    % Split up downstream elements as far as lDecay with logarithmic
    % spacing between elements
    splitInd=[]; splitBL={}; splitDS=[];
    ind1=lastEle+1;
    if ind1<=length(BEAMLINE)
      for ibl=ind1:length(BEAMLINE)
        if (BEAMLINE{ibl}.S-BEAMLINE{ind1}.S)>lDecay || strcmp(BEAMLINE{ibl}.Class,'SBEN')
          break
        end
        if isfield(BEAMLINE{ibl},'L') && BEAMLINE{ibl}.L>0
          splitInd(end+1)=ibl;
          splitBL{end+1}=BEAMLINE{ibl};
          splitDS(end+1)=BEAMLINE{ibl}.S-BEAMLINE{ind1}.S+BEAMLINE{ibl}.L;
        end
      end
    end
    if isempty(splitInd); iele=lastEle+1; continue; end;
    dl=logspace(-3,0,dsplit).*lDecay;
%     dl=linspace(0,lDecay,dsplit);
    blGrowInd=0;
    for isele=1:length(splitInd)
      BL=splitBL{isele}; olen=BL.L; BLorig=BL;
      if isele==1
        ndl=find(dl<=splitDS(isele));
      else
        ndl=find(dl>splitDS(isele-1) & dl<=splitDS(isele));
      end
      if isempty(ndl) || length(ndl)<2; continue; end;
      newBL={};
      for idl=1:length(ndl)
        if idl==1 && isele==1
          BL.L=dl(ndl(idl));
        elseif idl==1
          BL.L=dl(ndl(idl))-splitDS(isele-1);
        elseif idl==length(ndl)
          BL.L=splitDS(isele)-dl(ndl(idl-1));
        else
          BL.L=dl(ndl(idl))-dl(ndl(idl)-1);
        end
        BL.TrackFlag.doCSR=true;
        if isfield(BL,'B')
          BL.B=BLorig.B.*(BL.L/olen);
        end
        newBL{end+1}=BL;
      end
      if (splitInd(isele)+blGrowInd)<length(BEAMLINE)
        BEAMLINE=[BEAMLINE(1:splitInd(isele)+blGrowInd-1); newBL'; BEAMLINE(splitInd(isele)+1+blGrowInd:end)];
      else
        BEAMLINE=[BEAMLINE(1:splitInd(isele)+blGrowInd-1); newBL'];
      end
      blGrowInd=blGrowInd+length(newBL)-1;
    end
    lastEle=splitInd(end)+blGrowInd;
    SetSPositions( 1, length(BEAMLINE), 0 );
  end
  iele=lastEle+1;
end