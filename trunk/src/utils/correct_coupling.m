function [stat Bcorrect skewps] = correct_coupling(otruse, BEAMLINE, PS, GIRDER, FL)
stat{1}=1;


%intensities to scan to obtain response matrix
intensities=[-18 -10 -5 0 5 10 18];

%intensities to B
IBlookup=[-20 -16 -12 -8 -4 0 4 8 12 16 20;
         -0.2227 -0.1779 -0.1331 -0.0884 -0.0439 0 0.0439 0.0884 0.1331 0.1779 0.2227];
B = interp1(IBlookup(1,:),IBlookup(2,:),intensities);

%Find OTRs
for i = 0:3
  otrind(i+1) = findcells(BEAMLINE,'Name',['OTR',num2str(i),'X']);
end

%Find skew indices
for i = 1:4
  temp = findcells(BEAMLINE,'Name',['QK',num2str(i),'X']);
  skewind(i)=temp(2); 
end

%find the skewpowersupplies indexs
 for n=1:4
  skewps(n)=BEAMLINE{skewind(n)}.PS;
 end
 
 
%-------------------------
%GET RESPONSE MATRIX R
%------------------------

latticeversion=FL.SimModel.opticsVersion;
latticename=FL.SimModel.opticsName;
filename=sprintf('RespMatrix4CoupCorr_%s_%s',latticeversion,latticename);

%load response matrix it if exists
if exist(sprintf('/home/atf2-fs/ATF2/FlightSim/userData/%s.mat',filename),'file')
  load(sprintf('/home/atf2-fs/ATF2/FlightSim/userData/%s.mat',filename),'R');
else

Beam1=MakeBeam6DGauss(FL.SimModel.Initial,10001,5,0);

%Track the beam from element 1 to IEX
[stat,beam_EXT]=TrackThru(1,FL.SimModel.extStart,Beam1,1,1,0);

%GET for different B for each of the 4 skews the coupling terms in each OTR
for i=1:4%for the different skews
    for j=1:length(B)%for the different B
        
        %all skews off except i
        for n=1:4
            PS(skewps(n)).Ampl=0;
            PS(skewps(n)).SetPt=0;
        end
        PS(skewps(i)).Ampl=1;
        PS(skewps(i)).SetPt=1;

        %put B in the skew
        BEAMLINE{skewind(i)-1}.B=B(j);
        BEAMLINE{skewind(i)}.B=B(j); 

        %track the beam from EXT to all 4 OTRs
        [stat,beam_OTR0]=TrackThru(FL.SimModel.extStart,otrind(1),beam_EXT,1,1,0);
        [x,sigma_OTR0] = GetBeamPars(beam_OTR0,1);sigma_OTR0=sigma_OTR0(1:4,1:4);
        [stat,beam_OTR1]=TrackThru(otrind(1),otrind(2),beam_OTR0,1,1,0);
        [x,sigma_OTR1] = GetBeamPars(beam_OTR1,1);sigma_OTR1=sigma_OTR1(1:4,1:4);
        [stat,beam_OTR2]=TrackThru(otrind(2),otrind(3),beam_OTR1,1,1,0);
        [x,sigma_OTR2] = GetBeamPars(beam_OTR2,1);sigma_OTR2=sigma_OTR2(1:4,1:4);
        [stat,beam_OTR3]=TrackThru(otrind(3),otrind(4),beam_OTR2,1,1,0);
        [x,sigma_OTR3] = GetBeamPars(beam_OTR3,1);sigma_OTR3=sigma_OTR3(1:4,1:4);

        coup(i,1,j)=sigma_OTR0(1,3);
        coup(i,2,j)=sigma_OTR1(1,3);
        coup(i,3,j)=sigma_OTR2(1,3);
        coup(i,4,j)=sigma_OTR3(1,3);
       
    end
end

 %Do linear fit and build response matrix
 for i=1:4 %for each skew quad
     for k=1:4
       temp(1:length(B))=coup(i,k,:);
       c=polyfit(B(:),temp(:),1);
       R(k,i)=c(1);
     end    
 end
save(sprintf('/home/atf2-fs/ATF2/FlightSim/userData/%s.mat',filename),'R')
end
 
 l=0;
 for i=find(otruse)
    l=l+1;
    [stat,data]=getOTRsize(i);
    if (stat{1}~=1)
      stat{1}=-1;
      stat{2}=sprintf('Failed to get data for OTR%dX',n-1);
      return
    end
    sig13(l)=data.sig13;
 end

 sig13=sig13*1e-12;%stored data was in um
  
 if any(isnan(sig13))
      stat{1}=-1;
      stat{2}=sprintf('getOTRsize gave a NaN return');
      return
 end
 
%Calculate INT to correct
Bcorrect=R(find(otruse),:)\sig13';

%if any(abs(Bcorrect)>0.2227)%if intensity exceeds +-20A put +-20A
%   Bcorrect(Bcorrect>0.2227)=0.2227;
%   Bcorrect(Bcorrect<-0.2227)=-0.2227;
%end
