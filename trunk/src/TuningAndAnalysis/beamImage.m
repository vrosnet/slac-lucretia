function data=beamImage(beam,nsig,E0,asym,nbins,axhan,dpk,numpix)
% data=beamImage(beam,nsig [,E0,asym,nbins,axhan,dpk])
%
% Graphical plot of Lucretia beam in transverse plane including gaussian fits
% beam = Lucretia beam structure
% nsig = cut on transverse dimension
% E0 (optional) = centroid beam energy
% asym (optional) = true (default) means asymmetric fit
% nbins (optional) = number of histogram bins to use (default 150)
% axhan (optional) = supply axes handles to plot [trans long], if ==0 then
%                  omit
% dpk (optional) = perform double-peak fit for [trans long]
% numpix (optional) = number of pixels to plot for 2D plots [default using
%                     nbins^2]
% ===========================================
% GW, Sept 15, 2013
%   - add longitudinal and options
% MDW, Oct 10 2012

% Can we use GPU?
try
  if gpuDeviceCount>0
    dogpu=true;
  else
    dogpu=false;
  end
catch
  dogpu=false;
end

% Turn off frequent annoying warnings
w1=warning('query','MATLAB:rankDeficientMatrix');
w2=warning('query','MATLAB:nearlySingularMatrix');
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix

if ~exist('dpk','var') || isempty(dpk); dpk=[false false]; end;
if ~exist('axhan','var'); axhan=[]; end;

id=find(beam.Bunch.stop==0);
rays=beam.Bunch.x(:,id)'; %#ok<*FNDSB>
Q=beam.Bunch.Q(id);

conv=[1e6,1e6,1e6,1e6,1e6,1e2];
if ~exist('nbins','var') || isempty(nbins); nbins=150; end;
if ~exist('asym','var') || isempty(asym); asym=false; end;
if ~exist('numpix','var') || isempty(numpix)
  npixsqr=nbins;
else
  npixsqr=floor(sqrt(numpix));
end

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
  for n=1:length(v)
    eval(['u=',v{n},';'])
    N=length(u);
    nit=0;
    Qtmp=Q;
    while 1
      u0=mean(u);
      sig=sqrt(var(u,Qtmp));
      id=find(abs(u-u0)<nsig*sig);
      u=u(id);
      Qtmp=Q(id);
      if (length(u)==N),break,end
      N=length(u);
      nit=nit+1;
    end
    eval(sprintf('id%s=find(ismember(%s,u));',v{n},v{n}));
  end
else
  idx=(1:length(x))';
  idy=(1:length(x))';
  idz=(1:length(x))';
  iddp=(1:length(x))';
  idpx=(1:length(x))'; %#ok<NASGU>
  idpy=(1:length(x))'; %#ok<NASGU>
end


if isempty(axhan); axhan(1)=figure; axhan(2)=figure; end;

if axhan(1)~=0
  sh=subplot(2,2,1,'Parent',axhan(1));
  u=linspace(min(y(idy)),max(y(idy)),nbins);
  [~,BIN]=histc(y(idy),u); v=accumarray(BIN,Q(idy)).*1e9;
  bar(sh,u,v);
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
  data.sigy=q(4);
  h1=barh(sh,u,v);
  set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
  hold(sh,'on')
  plot(sh,yfit,u,'r-')
  hold(sh,'off')
  title(sh,sprintf('\\sigma_y = %.1f um',abs(data.sigy)))
  ylabel(sh,'y (um)')
  xlabel(sh,'Q (nC)')
  sh_y=sh;

  sh=subplot(2,2,4,'Parent',axhan(1));
  u=linspace(min(x(idx)),max(x(idx)),nbins);
  [~,BIN]=histc(x(idx),u); v=accumarray(BIN,Q(idx)).*1e9;
  bar(sh,u,v);
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
  data.sigx=q(4);
  h2=bar(sh,u,v);
  set(h2,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
  hold(sh,'on')
  plot(sh,u,xfit,'r-')
  hold(sh,'off')
  title(sh,sprintf('\\sigma_x = %.1f um',abs(data.sigx)))
  xlabel(sh,'x (um)')
  ylabel(sh,'Q (nC)')
  sh_x=sh;

  sh=subplot(2,2,1,'Parent',axhan(1));axis(sh,'square');
  sh=subplot(2,2,4,'Parent',axhan(1));axis(sh,'square');

  sh=subplot(2,2,2,'Parent',axhan(1));
  if dogpu
    id=idx(ismember(gpuArray(idx),gpuArray(idy)));
  else
    id=idx(ismember(idx,idy));
  end
  data.rmsx=std(x(idx));
  data.rmsy=std(y(idy));
  xb=linspace(min(x(id)),max(x(id)),npixsqr);
  yb=linspace(min(y(id)),max(y(id)),npixsqr);
  xr=interp1(xb,1:npixsqr,x(id),'nearest');
  yr=interp1(yb,1:npixsqr,y(id),'nearest');
  beamprof = accumarray([xr(:) yr(:)], Q(id), [npixsqr npixsqr]);
  imagesc(xb,yb,beamprof','Parent',sh);
  set(sh,'YDir','normal')
  axis(sh,'square');
  ax=axis(sh);
  ax_x=axis(sh_x); ax_x(1:2)=ax(1:2); axis(sh_x,ax_x);
  ax_y=axis(sh_y); ax_y(3:4)=ax(3:4); axis(sh_y,ax_y);
  ylabel(sh,'y (um)')
  xlabel(sh,'x (um)')
  
  sh=subplot(2,2,3,'Parent',axhan(1));
  cla(sh)
  axis(sh,'off');
  text(0.1,0.8,sprintf('rms X = %g',data.rmsx),'Parent',sh);
  text(0.1,0.6,sprintf('rms Y = %g',data.rmsy),'Parent',sh);
  text(0.1,0.4,sprintf('Q = %g nC',1e9*sum(beam.Bunch.Q(~beam.Bunch.stop))),'Parent',sh)
end

if length(axhan)==1 || axhan(2)==0; return; end;

sh=subplot(2,2,1,'Parent',axhan(2));
u=linspace(min(dp(iddp)),max(dp(iddp)),nbins);
[~,pbins]=histc(dp(iddp),u); v=accumarray(pbins,Q(iddp)).*1e9;
bar(sh,u,v);
if (asym)
  [yfit,q]=agauss_fit(u,v,[],0);
  xi=u;
else
  if dpk(2)
    [q,~,~,xi,yfit]=peakfit([u' v],0,0,2,1,0,10,0,0,0,0); qE=q;
  else
    [yfit,q]=gauss_fit(u,v,[],0);
    xi=u;
  end
end
data.sigE=q(4);
h1=barh(sh,u,v);
set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold(sh,'on')
plot(sh,yfit,xi,'r-')
hold(sh,'off')
if ~dpk(2); title(sh,sprintf('\\sigma_{dP} = %.1g %%',abs(data.sigE))); end;
xlabel(sh,'Q (nC)')
ylabel(sh,'dP (%)')
sh_p=sh;

sh=subplot(2,2,4,'Parent',axhan(2));
cla(sh)
u=linspace(min(z(idz)),max(z(idz)),nbins);
[~,zbins]=histc(z(idz),u); v=accumarray(zbins,Q(idz));
v=1e-3.*v.*(1/(((u(2)-u(1))*1e-6)/299792458)); % y-axis Q->I (kA)
bar(sh,u,v);
if (asym)
  [xfit,q]=agauss_fit(u,v,[],0);
  xi=u;
else
  if dpk(2)
    [q,~,~,xi,xfit]=peakfit([u v'],0,0,2,1,0,10,0,0,0,0);
  else
    [xfit,q]=gauss_fit(u,v,[],0);
    xi=u;
  end
end
data.sigz=q(4);
h2=bar(sh,u,v);
set(h2,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold(sh,'on')
plot(sh,xi,xfit,'r-')
hold(sh,'off');
if ~dpk(2); title(sh,sprintf('\\sigma_z = %.1f um I(pk) = %.1f kA',abs(data.sigz),max(v))); end;
data.pkI=max(v);
xlabel(sh,'z (um)')
ylabel(sh,'I (kA)');
sh_z=sh;

sh=subplot(2,2,1,'Parent',axhan(2));axis(sh,'square');
sh=subplot(2,2,4,'Parent',axhan(2));axis(sh,'square');

sh=subplot(2,2,2,'Parent',axhan(2));
if dogpu
  id=idz(ismember(gpuArray(idz),gpuArray(iddp)));
else
  id=idz(ismember(idz,iddp));
end
data.rmsz=std(z(idz));
data.rmsE=std(dp(iddp));
zb=linspace(min(z(id)),max(z(id)),npixsqr);
pb=linspace(min(dp(id)),max(dp(id)),npixsqr);
zr=interp1(zb,1:npixsqr,z(id),'nearest');
pr=interp1(pb,1:npixsqr,dp(id),'nearest');
beamprof = accumarray([zr(:) pr(:)], Q(id), [npixsqr npixsqr]);
imagesc(zb,pb,beamprof','Parent',sh);

sh=subplot(2,2,2,'Parent',axhan(2));
set(sh,'YDir','normal')
axis(sh,'square');
ax=axis(sh);
ax_x=axis(sh_z); ax_x(1:2)=ax(1:2); axis(sh_z,ax_x);
ax_y=axis(sh_p); ax_y(3:4)=ax(3:4); axis(sh_p,ax_y);
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
  text(0,0.8,sprintf('rms Z = %g',data.rmsz),'Parent',sh);
end

drawnow('expose')

warning(w1.state,'MATLAB:rankDeficientMatrix');
warning(w2.state,'MATLAB:nearlySingularMatrix');
