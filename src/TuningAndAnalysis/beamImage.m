function [sigx,sigy,rmsx,rmsy,pkI]=beamImage(beam,nsig,E0,alim)
% [sigx,sigy,rmsx,rmsy]=beamImage(beam,nsig [,E0,alim])
%
% Graphical plot of Lucretia beam in transverse plane including gaussian fits
% beam = Lucretia beam structure
% nsig = cut on transverse dimension
% E0 (optional) = centroid beam energy
% alim (optional) = force axis limits
%
% ===========================================
% GW, Sept 15, 2013
%   - add longitudinal and options
% MDW, Oct 10 2012

id=find(beam.Bunch.stop==0);
rays=beam.Bunch.x(:,id)';

name={'X (um)','PX (ur)','Y (um)','PY (ur)','Z (um)','dP (%)'};
conv=[1e6,1e6,1e6,1e6,1e6,1e2];
nbin=150;
asym=true;
npix=250;

x=conv(1)*rays(:,1);
px=conv(2)*rays(:,2);
y=conv(3)*rays(:,3);
py=conv(4)*rays(:,4);
z=conv(5)*rays(:,5);
E=rays(:,6);
if (~exist('E0','var')),E0=mean(E);end
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
  idpx=(1:length(x))';
  idy=(1:length(x))';
  idpy=(1:length(x))';
  idz=(1:length(x))';
end
iddp=(1:length(x))';

figure
subplot(221)
[v,u]=hist(y,nbin);
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix
if (asym)
  [yfit,q,dq,chi2]=agauss_fit(u,v,[],0);
else
  [yfit,q,dq,chi2]=gauss_fit(u,v,[],0);
end
warning on MATLAB:rankDeficientMatrix
warning on MATLAB:nearlySingularMatrix
sigy=q(4);
h1=barh(u,v);
ylim=get(gca,'YLim');
set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold on
plot(yfit,u,'r-')
hold off
title(sprintf('\\sigma_y = %.1f um',abs(sigy)))
ylabel('y (um)')

subplot(224)
[v,u]=hist(x,nbin);
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix
if (asym)
  [xfit,q,dq,chi2]=agauss_fit(u,v,[],0);
else
  [xfit,q,dq,chi2]=gauss_fit(u,v,[],0);
end
warning on MATLAB:rankDeficientMatrix
warning on MATLAB:nearlySingularMatrix
sigx=q(4);
h2=bar(u,v);
xlim=get(gca,'XLim');
set(h2,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold on
plot(u,xfit,'r-')
hold off
title(sprintf('\\sigma_x = %.1f um',abs(sigx)))
xlabel('x (um)')

ulim=max(abs([xlim,ylim]));
if (~exist('alim','var'))
  alim=ulim;
end
subplot(221),set(gca,'YLim',[-alim,alim]),axis square
subplot(224),set(gca,'XLim',[-alim,alim]),axis square

subplot(222)
id=intersect(idx,idy);
rmsx=std(x(id));
rmsy=std(y(id));
beamprof=hist2(x(id),y(id), ...
  linspace(-ulim,ulim,npix),linspace(-ulim,ulim,npix));
imagesc([-ulim,ulim],[-ulim,ulim],beamprof)
set(gca,'YDir','normal')
axis square
axis([-alim,alim,-alim,alim])
title(sprintf('%d\\sigma cuts: rmsx= %.1f um, rmsy= %.1f um', ...
  nsig,rmsx,rmsy))
ylabel('y (um)')
xlabel('x (um)')

figure
subplot(221)
[v,u]=hist(dp,nbin);
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix
if (asym)
  [yfit,q,dq,chi2]=agauss_fit(u,v,[],0);
else
  [yfit,q,dq,chi2]=gauss_fit(u,v,[],0);
end
warning on MATLAB:rankDeficientMatrix
warning on MATLAB:nearlySingularMatrix
sigE=q(4);
h1=barh(u,v);
ylim=get(gca,'YLim');
set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold on
plot(yfit,u,'r-')
hold off
title(sprintf('\\sigma_E = %.1f %% RMS = %.1f %%',abs(sigE),100*(std(E)/mean(E))))
ylabel('dP (%)')

subplot(224)
% [v,u]=hist(z(z>337),nbin);
[v,u]=hist(z,nbin);
dQ=sum(beam.Bunch.Q(~beam.Bunch.stop))/sum(~beam.Bunch.stop);
v=1e-3.*v.*(dQ/(((u(2)-u(1))*1e-6)/299792458));
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix
if (asym)
  [xfit,q,dq,chi2]=agauss_fit(u,v,[],0);
else
  [xfit,q,dq,chi2]=gauss_fit(u,v,[],0);
end
warning on MATLAB:rankDeficientMatrix
warning on MATLAB:nearlySingularMatrix
sigx=q(4);
h2=bar(u,v);
xlim=get(gca,'XLim');
set(h2,'EdgeColor',[0,0,1],'FaceColor',[0,0,1])
hold on
plot(u,xfit,'r-')
hold off
title(sprintf('\\sigma_z = %.1f um I(pk) = %.1f kA',abs(sigx),max(v)))
pkI=max(v);
xlabel('z (um)')
ylabel('I (kA)');

ulim=max(abs([xlim,ylim]));
if (~exist('alim','var'))
  alim=ulim;
end
% subplot(221),set(gca,'YLim',[-alim,alim]),axis square
% subplot(224),set(gca,'XLim',[-alim,alim]),axis square

subplot(222)
id=intersect(idz,iddp);
beamprof=hist2(z(id),dp(id), ...
  linspace(xlim(1),xlim(2),npix),linspace(ylim(1),ylim(2),npix));
imagesc(xlim,ylim,beamprof)
% plot(z(id),dp(id),'.')
set(gca,'YDir','normal')
% axis square
axis([xlim ylim])
title(sprintf('Mean Energy = %.3f GeV',mean(E)))
ylabel('dP (%)')
xlabel('z (um)')

drawnow('expose')

end
