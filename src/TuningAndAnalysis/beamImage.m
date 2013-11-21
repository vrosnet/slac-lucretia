function [sigx,sigy,rmsx,rmsy,pkI]=beamImage(beam,nsig,E0,asym,nbins,axhan,dpk,numpix)
% [sigx,sigy,rmsx,rmsy]=beamImage(beam,nsig [,E0,asym,nbins,axhan,dpk])
%
% Graphical plot of Lucretia beam in transverse plane including gaussian fits
% beam = Lucretia beam structure
% nsig = cut on transverse dimension
% E0 (optional) = centroid beam energy
% asym (optional) = true (default) means asymmetric fit
% nbins (optional) = number of histogram bins to use (default 200)
% axhan (optional) = supply axes handles to plot [trans long], if ==0 then
%                  omit
% dpk (optional) = perform double-peak fit for [trans long]
% numpix (optional) = number of pixels to plot for 2D plots
% ===========================================
% GW, Sept 15, 2013
%   - add longitudinal and options
% MDW, Oct 10 2012
w1=warning('query','MATLAB:rankDeficientMatrix');
w2=warning('query','MATLAB:nearlySingularMatrix');
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix

if ~exist('dpk','var') || isempty(dpk); dpk=[false false]; end;
if ~exist('axhan','var'); axhan=[]; end;

id=find(beam.Bunch.stop==0);
rays=beam.Bunch.x(:,id)'; %#ok<*FNDSB>

conv=[1e6,1e6,1e6,1e6,1e6,1e2];
if ~exist('nbins','var') || isempty(nbins); nbins=200; end;
nbin=nbins;
if ~exist('asym','var') || isempty(asym); asym=false; end;
if ~exist('numpix','var') || isempty(numpix); numpix=250^2; end;
npix=floor(sqrt(numpix));

x=conv(1)*rays(:,1);
px=conv(2)*rays(:,2); %#ok<NASGU>
y=conv(3)*rays(:,3);
py=conv(4)*rays(:,4); %#ok<NASGU>
z=conv(5)*rays(:,5);
E=rays(:,6);
if (~exist('E0','var')) || isempty(E0); E0=mean(E);end
dp=conv(6)*(E-E0)/E0;

if (nsig>0)
  v={'x','px','y','py','z','dp'};
  for n=1:5 % don't do dp (FACET fully compressed)
    eval(['u=',v{n},';'])
    N=length(u);
    nit=0;
    while 1
      u0=mean(u);
      sig=std(u);
      id=find(abs(u-u0)<nsig*sig);
      u=u(id);
      if (length(u)==N),break,end
      N=length(u);
      nit=nit+1;
    end
    eval([v{n},'=u(id);'])
    eval(['id',v{n},'=id;'])
  end
else
  idx=(1:length(x))';
  idy=(1:length(x))';
  idz=(1:length(x))';
end
iddp=(1:length(x))';

if isempty(axhan); axhan(1)=figure; axhan(2)=figure; end;

if axhan(1)~=0
  sh=subplot(2,2,1,'Parent',axhan(1));
  [v,u]=hist(sh,y,nbin);
  if (asym)
    if dpk(1)
      [yfit,q]=agauss_fit2(u,v,[],0);
    else
      [yfit,q]=agauss_fit(u,v,[],0);
    end
  else
    if dpk(1)
      [yfit,q]=gauss_fit(u,v,[],0);
    else
      [yfit,q]=gauss_fit(u,v,[],0);
    end
  end
  sigy=q(4);
  h1=barh(sh,u,v);
  ylim=get(sh,'YLim');
  set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
  hold(sh,'on')
  plot(sh,yfit,u,'r-')
  hold(sh,'off')
  title(sh,sprintf('\\sigma_y = %.1f um',abs(sigy)))
  ylabel(sh,'y (um)')

  sh=subplot(2,2,4,'Parent',axhan(1));
  [v,u]=hist(sh,x,nbin);
  if (asym)
    if dpk(1)
      [xfit,q]=agauss_fit2(u,v,[],0);
    else
      [xfit,q]=agauss_fit(u,v,[],0);
    end
  else
    if dpk(1)
      [xfit,q]=gauss_fit2(u,v,[],0);
    else
      [xfit,q]=gauss_fit(u,v,[],0);
    end
  end
  sigx=q(4);
  h2=bar(sh,u,v);
  xlim=get(sh,'XLim');
  set(h2,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
  hold(sh,'on')
  plot(sh,u,xfit,'r-')
  hold(sh,'off')
  title(sh,sprintf('\\sigma_x = %.1f um',abs(sigx)))
  xlabel(sh,'x (um)')

  ulim=max(abs([xlim,ylim]));
  sh=subplot(2,2,1,'Parent',axhan(1));axis(sh,'square');
  sh=subplot(2,2,4,'Parent',axhan(1));axis(sh,'square');

  sh=subplot(2,2,2,'Parent',axhan(1));
  id=intersect(idx,idy);
  rmsx=std(x(id));
  rmsy=std(y(id));
  beamprof=hist2(x(id),y(id), ...
    linspace(-ulim,ulim,npix),linspace(-ulim,ulim,npix));
  imagesc([-ulim,ulim],[-ulim,ulim],beamprof,'Parent',sh)
  set(sh,'YDir','normal')
  axis(sh,'square');
  ylabel(sh,'y (um)')
  xlabel(sh,'x (um)')
  
  sh=subplot(2,2,3,'Parent',axhan(1));
  cla(sh)
  axis(sh,'off');
  text(0.1,0.8,sprintf('rms X = %g',rmsx),'Parent',sh);
  text(0.1,0.6,sprintf('rms Y = %g',rmsy),'Parent',sh);
  text(0.1,0.4,sprintf('Q = %g nC',1e9*sum(beam.Bunch.Q(~beam.Bunch.stop))),'Parent',sh)
end

if length(axhan)==1 || axhan(2)==0; return; end;

sh=subplot(2,2,1,'Parent',axhan(2));
[v,u]=hist(sh,dp,nbin);
if (asym)
  [yfit,q]=agauss_fit(u,v,[],0);
  xi=u;
else
  if dpk(2)
    [q,~,~,xi,yfit]=peakfit([u' v'],0,0,2,1,0,10,0,0,0,0); qE=q;
  else
    [yfit,q]=gauss_fit(u,v,[],0);
    xi=u;
  end
end
sigE=q(4);
h1=barh(sh,u,v);
ylim=get(sh,'YLim');
set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold(sh,'on')
plot(sh,yfit,xi,'r-')
hold(sh,'off')
% title(sh,sprintf('\\sigma_E = %.1f %% RMS = %.1f %%',abs(sigE),100*(std(E)/mean(E))))
ylabel(sh,'dP (%)')

sh=subplot(2,2,4,'Parent',axhan(2));
cla(sh)
[v,u]=hist(sh,z,nbin);
dQ=sum(beam.Bunch.Q(~beam.Bunch.stop))/sum(~beam.Bunch.stop);
v=1e-3.*v.*(dQ/(((u(2)-u(1))*1e-6)/299792458));
if (asym)
  [xfit,q]=agauss_fit(u,v,[],0);
  xi=u;
else
  if dpk(2)
    [q,~,~,xi,xfit]=peakfit([u' v'],0,0,2,1,0,10,0,0,0,0);
  else
    [xfit,q]=gauss_fit(u,v,[],0);
    xi=u;
  end
end
sigz=q(4);
h2=bar(sh,u,v);
xlim=get(sh,'XLim');
set(h2,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold(sh,'on')
plot(sh,xi,xfit,'r-')
hold(sh,'off');
if ~dpk(2); title(sh,sprintf('\\sigma_z = %.1f um I(pk) = %.1f kA',abs(sigz),max(v))); end;
pkI=max(v);
xlabel(sh,'z (um)')
ylabel(sh,'I (kA)');

sh=subplot(2,2,1,'Parent',axhan(2));axis(sh,'square');
sh=subplot(2,2,4,'Parent',axhan(2));axis(sh,'square');

sh=subplot(2,2,2,'Parent',axhan(2));
id=intersect(idz,iddp);
rmsz=std(z(id));
beamprof=hist2(z(id),dp(id), ...
  linspace(xlim(1),xlim(2),npix),linspace(ylim(1),ylim(2),npix));
imagesc(xlim,ylim,beamprof,'Parent',sh)
sh=subplot(2,2,2,'Parent',axhan(2));axis(sh,'square');
set(sh,'YDir','normal')
axis(sh,[xlim ylim])
title(sh,sprintf('Mean Energy = %.3f GeV',mean(E)))
ylabel(sh,'dP (%)')
xlabel(sh,'z (um)')

sh=subplot(2,2,3,'Parent',axhan(2));
cla(sh)
axis(sh,'off');

if dpk(2)
  text(-0.2,0.8,sprintf('pk 1 (pos, ampl, width):'),'Parent',sh);
  text(-0.2,0.7,sprintf('[Z]: %.3g, %.3g, %.3g',q(1,2),q(1,3),q(1,4)),'Parent',sh);
  text(-0.2,0.6,sprintf('[E]: %.3g, %.3g, %.3g',qE(1,2),qE(1,3),qE(1,4)),'Parent',sh);
  text(-0.2,0.5,sprintf('pk 2 (pos, ampl, width):'),'Parent',sh);
  text(-0.2,0.4,sprintf('[Z]: %.3g, %.3g, %.3g',q(2,2),q(2,3),q(2,4)),'Parent',sh);
  text(-0.2,0.3,sprintf('[E]: %.3g, %.3g, %.3g',qE(2,2),qE(2,3),qE(2,4)),'Parent',sh);
  text(-0.2,0.1,sprintf('dZ = %.3g um',abs(q(1,2)-q(2,2))),'Parent',sh);
else
  text(0,0.9,sprintf('rms E = %g',100*(std(E)/mean(E))),'Parent',sh);
  text(0,0.8,sprintf('rms Z = %g',rmsz),'Parent',sh);
end

drawnow('expose')

warning(w1.state,'MATLAB:rankDeficientMatrix');
warning(w2.state,'MATLAB:nearlySingularMatrix');
