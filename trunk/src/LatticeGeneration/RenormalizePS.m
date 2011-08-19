function stat = RenormalizePS( klysno )

% RENORMALIZEPS Rescale a power supply such that its amplitude == 1.0
%
%    stat = RenormalizePS( PSno ) rescales the Ampl, SetPt, and Step of a
%       power supply such that Ampl is returned to 1.0; the inverse of the
%       amplitude scale factor is applied to each magnet assigned to the
%       PS.  Return variable stat is a cell array, with stat{1} == 1 if the
%       exchange occurred without error, == 0 if errors occurred, and
%       stat{2...} are text error messages.

% MOD:
%        29-sep-2005, PT:
%           support for magnets with multiple power supplies.
%        28-sep-2005, PT:
%           improved handling of zero amplitude.

%==========================================================================

  global BEAMLINE ;
  global PS ;
  stat = InitializeMessageStack( ) ;
  
% Is the desired PS in range?

  if ( klysno > length(PS) )
      stat = AddMessageToStack(stat,...
          ['PS # ',num2str(klysno),...
          ' out of range in MovePhysicsVarsToPS']) ;
      stat{1} = 0 ;
  end
  
% if the PS amplitude is zero, just go to its magnets and set their B
% values to 1

  if (PS(klysno).Ampl == 0)
     for elemno = PS(klysno).Element
       if (length(BEAMLINE{elemno}.PS) == 1)
         BEAMLINE{elemno}.B = ones(1,length(BEAMLINE{elemno}.B)) ;
       else   
         PS_index = find(BEAMLINE{elemno}.PS == klysno) ;
         BEAMLINE{elemno}.B(PS_index) = 1 ;
       end    
     end
     return ;
  end
  
% compute the scale factor

  scale = 1 / PS(klysno).Ampl ;
  
% apply the scale factor to the PS

  PS(klysno).Ampl = PS(klysno).Ampl * scale ;
  PS(klysno).Step = PS(klysno).Step * scale ;
  PS(klysno).SetPt = PS(klysno).SetPt * scale ;
  
% now apply the reverse transformation on elements 

  for elemno = PS(klysno).Element
       if (length(BEAMLINE{elemno}.PS) == 1) 
         BEAMLINE{elemno}.B = BEAMLINE{elemno}.B / scale ;
       else
         PS_index = find(BEAMLINE{elemno}.PS == klysno) ;  
         BEAMLINE{elemno}.B(PS_index) = BEAMLINE{elemno}.B(PS_index) / scale ;
       end    
  end