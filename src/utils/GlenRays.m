% load 'beam.mat'
rays0=[beamout.Bunch.x]';
P=mean(rays0(:,6)); % GeV/c
dp=(rays0(:,6)-P)/P;
rays0(:,6)=dp;

[rays,xv,sigm,sigv,emit,emitn,twss,eta1,eta2]=AnalyzeRays(rays0,P,3,1);
PlotRays(rays,'Lucretia',0)
