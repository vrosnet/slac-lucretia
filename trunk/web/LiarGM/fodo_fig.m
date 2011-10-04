
sbpm1m=sbpm1-2560.;
sbpm2m=2560.-sbpm2;

for i=1:64
    strr=num2str(0.01*(i-1))
    timestr= ['time (s) = ',strr];
    figure
  hs=plot(sbpm1m,mislx1(:,i),'b')
    hold on
  hs=plot(sbpm2m,mislx2(:,i),'g')
  hs=plot(sbpm1m,beamx1(:,i),'r')
  hs=plot(sbpm2m,beamx2(:,i),'m')  
    
  hc=title('Two FODO linacs')
    set(hc,'FontSize',14)  
  hc=ylabel('x_{ground}, x_{beam}, m')
    set(hc,'FontSize',12)
  hc=xlabel('s, m')
    set(hc,'FontSize',12)  
  axis([-3000 3000 -6.e-8 6.e-8 ]) 
  text(1500, 5.e-8, timestr)
  legend('x1_{ground}','x2_{ground}','x1_{beam}','x2_{beam}',2)
  pause
  close
end 


figure
hc=plot(beamx1(128,:)-beamx2(128,:),'.')
hc=title('e+ e- beam separation. FODO example')
set(hc,'FontSize',14)
hc=ylabel('beam separation, m')
set(hc,'FontSize',12)
hc=xlabel('pulse number (@100Hz)')
set(hc,'FontSize',12)
hc=axis([0 70 -6.e-8 6.e-8 ])

