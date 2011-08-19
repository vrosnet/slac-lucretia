function stat = PSTrim( klyslist, varargin )

% PSTRIM Set power supply actual values to desired values
%
%   stat = PSTrim( PS_List ) sets the Ampl field of each PS in the PS_List
%      equal to the SetPt of that device.  The step sizes are taken into
%      account if they are not zero, which can result in a PS with a
%      residual difference in its desired and actual strength. Return
%      argument stat is a cell array, with stat{1} == 1 indicating success
%      and stat{1} == 0 indicating error.  Error messages are transferred
%      in stat{2...}.
%
% See also GirderMoverTrim, KlystronTrim.

%==========================================================================

  global PS ;
  stat = InitializeMessageStack( ) ;
  if (max(klyslist) > length(PS))
      stat{1} = 0 ;
      stat = AddMessageToStack(stat,...
          'Out-of-range power supplies found in PSTrim') ;
  end
  
% loop over power supplies

  for count = 1:length(klyslist) 
      klysno = klyslist(count) ;
      if (PS(klysno).Step == 0)
          PS(klysno).Ampl = PS(klysno).SetPt ;
      else
          nstep = round( (PS(klysno).SetPt - ...
                          PS(klysno).Ampl        ) / ...
                          PS(klysno).Step            ) ;
          PS(klysno).Ampl = PS(klysno).Ampl + ...
              nstep * PS(klysno).Step ;
      end
  end
      