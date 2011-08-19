function PlotTwiss(istart,iend,Twiss,Title,functions)
%
% PLOTTWISS Generate a standard Twiss plot with a magnet display overhead
%
% PlotTwiss( start, end, Twiss, Title, select ) produces a single- or
%   dual-y-axis plot of the betatron and dispersion functions, and appends
%   a magnet display over the plot.  Scalar arguments start and end
%   determine the boundaries of the plot in the BEAMLINE cell array, and
%   Twiss is the output of a GetTwiss function which contains the necessary
%   data for plotting.  Title is a text string title for the plot.
%   Argument select is a 3-vector which indicates which functions are to be
%   plotted:  if select(1) ~= 0, betatron functions are plotted; if
%   select(2) ~= 0, horizontal dispersion is plotted; if select(3) ~= 0,
%   vertical dispersion is plotted.
%
% See also:  GetTwiss, AddMagnetPlot.
%
% Version date:  22-Feb-2008.

% MOD:
%      22-Feb-2008, PT:
%         minor adjustment to plotting bounds -- need to plot 
%         to the end of the last element, which is 1 slot further
%         in the Twiss data compared to the element number.

  S = Twiss.S(istart:iend+1) ;
  b1 = [] ; b2 = [] ; e1 = [] ; e2 = [] ;
  if (functions(1))
      b1 = Twiss.betax(istart:iend+1) ;
      b2 = Twiss.betay(istart:iend+1) ;
  end
  
  if (functions(2))
      e1 = Twiss.etax(istart:iend+1) ;
  end
  if (functions(3))
      e2 = Twiss.etay(istart:iend+1) ;
  end
  
  if ( (isempty(e1)) & (isempty(e2)) )
      PlotBetas(S,b1,b2,1) ;
      ax = [] ;
  elseif ( (isempty(b1)) & (isempty(b2)) )
      PlotBetas(S,e1,e2,2) ;
      ax = [] ;
  else
      [ax,h0,h1] = PlotBetaEta(S,b1,b2,e1,e2) ;
  end
  
  if (~isempty(Title))
    title(Title) ;
  end
  [h2,h3] = AddMagnetPlot(istart,iend) ;
  AddTwissLegends(b1,b2,e1,e2,h2,ax) ;
  
%==========================================================================

function PlotBetas(S,b1,b2,BetaOrEta)

  YLabel = GetTwissLabel(b1,b2,BetaOrEta) ;
  % figure ;
  if (~isempty(b1))
    plot(S,b1) ;
    if (~isempty(b2))
        hold on ;
    end
  end
  if (~isempty(b2))
      hold on 
      plot(S,b2,'r--') ;
  end
  hold off
  xlabel('S position [m]') ;
  ylabel(YLabel) ;
  
%==========================================================================

function l = GetTwissLabel(f1,f2,BetaOrEta) ;

  if (BetaOrEta == 1)
      l = '\beta_{' ;
  else
      l = '\eta_{' ;
  end
  if (~isempty(f1))
      l = [l,'x'] ;
      if (~isempty(f2))
          l = [l,','] ;
      end
  end
  if (~isempty(f2))
      l = [l,'y'] ;
  end
  l = [l,'} [m]'] ;

%==========================================================================  

function [ax,h0,h1] = PlotBetaEta(S,b1,b2,e1,e2) ;

  % figure ;
  Y1Label = GetTwissLabel(b1,b2,1) ;
  Y2Label = GetTwissLabel(e1,e2,2) ;
  XLabel = 'S Position [m]' ;
  
  plot(S,b1) ;
  hold on 
  plot(S,b2,'r--') ;
  xlabel(XLabel) ;
  ylabel(Y1Label) ;
  if (isempty(e1))
    [ax,h0,h1] = plotyy(S,b1,S,e2) ;
  else
    [ax,h0,h1] = plotyy(S,b1,S,e1) ;
  end
  if (  (~isempty(e1)) & (~isempty(e2))  ) 
    axes(ax(2)) ;
    hold on 
    plot(S,e2,'c--') ;
    hold off
  end
  ylabel(ax(2),Y2Label) ;
  
%==========================================================================  

function AddTwissLegends(b1,b2,e1,e2,h,ax)

  if (~isempty(b1))
      BetaLegend(h,ax) ;
  end
  if (  (~isempty(e1)) | (~isempty(e2))  )
      EtaLegend(e1,e2,h,ax) ;
  end
  
%==========================================================================  

function BetaLegend(h,ax)
  
  if (isempty(ax))
      axes(h)
      legend('\beta_x','\beta_y','Location','northwest')
  else
      legend(ax(1),'\beta_x','\beta_y','Location','northwest')
  end
  
%==========================================================================  

function EtaLegend(e1,e2,h,ax)

  if (isempty(ax))
      axes(h)
      if (isempty(e1))
        legend('\eta_y') ;
      elseif (isempty(e2))
         legend('\eta_x') ;
      else
         legend('\eta_x','\eta_y') ;
      end
  else
      if (isempty(e1))
        legend(ax(2),'\eta_y') ;
      elseif (isempty(e2))
         legend(ax(2),'\eta_x') ;
      else
         legend(ax(2),'\eta_x','\eta_y') ;
      end
  end
      
