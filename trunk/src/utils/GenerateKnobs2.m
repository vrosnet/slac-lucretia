function [x_knob, y_knob] = GenerateKnobs2(Group,model,beam0,Gind,PSind,Handles)

global BEAMLINE PS GIRDER VERB

v_init=VERB; VERB=0;

s_list=[1 2 4 5];
% for f=1:length(s_list)
%   s_bpmind(f)=Group.Sext.Bpm(s_list(f)).Han;
% end
k2=[ 1.69939443454 -0.419515957698 3.301438410805 -7.252873476491 7.740702674458];

% Get dispersion values at sextupoles
E_nominal=mean(beam0.Bunch.x(6,:));
disp_spread=2e-2;
ndisp=20;
disp=-disp_spread/2:disp_spread/(ndisp-1):disp_spread/2;
beam1=beam0;
for idisp=1:ndisp
  beam1.Bunch.x(6,:)=beam0.Bunch.x(6,:)+(E_nominal*disp(idisp));
  [stat,beamout,instdata] = TrackThru( 1, length(BEAMLINE), beam0, 1, 1, 0 );if stat{1}~=1; error(stat{2:end}); end;
  s_x(idisp,:)=[instdata{1}([Group.Sext.Bpm.Han]).x];
end
for ind=1:length(s_list);
  [q,dq,chisq,rms]=noplot_polyfit(disp,s_x(:,ind),1,1);
  eta_s(ind) = q(2);
end

for ind=1:length(s_list)
  isext=Group.Sext.dB.ClusterList(s_list(ind)).index(1);
  ds_x(ind)=k2(s_list(ind))*BEAMLINE{isext}.L*2*...
    model.Twiss.betax(isext)*model.Twiss.betax(model.ip_ind);
  ds_y(ind)=k2(s_list(ind))*BEAMLINE{isext}.L*2*...
    model.Twiss.betay(isext)*model.Twiss.betay(model.ip_ind);
  deta(ind)=k2(s_list(ind))*BEAMLINE{isext}.L*2*eta_s(ind)*...
    sqrt(model.Twiss.betax(isext)*model.Twiss.betax(model.ip_ind));
  M(:,ind)=[ds_x(ind); ds_y(ind); deta(ind);];
end

for f=1:3
  move=[0;0;0;];
  move(f)=1;
  x_knob{f}=lscov(M,move);
  x_knob{f}=x_knob{f}./max(abs(x_knob{f}));
end

s_init=GIRDER{Gind.Sext(5)}.MoverPos; bpm_ref=BEAMLINE{Group.Sext.Bpm(5).Ind}.ElecOffset;
ind=0;
smoves=[-1e-5:2e-5/9:1e-5];
for smove=smoves
  ind=ind+1;
  GIRDER{Gind.Sext(5)}.MoverSetPt(1)=s_init(2)+smove; stat=MoverTrim(Gind.Sext(5));
  BEAMLINE{Group.Sext.Bpm(5).Ind}.ElecOffset=bpm_ref+[0 smove];
%   [bs_mean, bs_std, errorlim] = BDS_steer(Group, PSind, beam0, model, Handles, Gind, 0.1);
  IPdat = IP_meas( Handles, PSind, beam0, model, [1 1] );
  s_disp(ind)=IPdat.Disp_y;
  s_coup(ind)=IPdat.sigma(2,3);
end
GIRDER{Gind.Sext(5)}.MoverSetPt=s_init; stat=MoverTrim(Gind.Sext(5)); BEAMLINE{Group.Sext.Bpm(5).Ind}.ElecOffset=bpm_ref;
[q,dq,chisq,rms]=noplot_polyfit(smoves,s_disp,1,1); ddispy_sm=q(2);
[q,dq,chisq,rms]=noplot_polyfit(smoves,s_coup,1,1); dcoup_sm=q(2);
p_init=PS(PSind.Quad(Group.SQuad_id(5))).Ampl;
ind=0;
pmoves=[-0.5:1/9:0.5];
for pmove=pmoves
  ind=ind+1;
  PS(PSind.Quad(Group.SQuad_id(5))).SetPt=p_init+pmove; stat=PSTrim(PSind.Quad(Group.SQuad_id(5)));
  [bs_mean, bs_std, errorlim] = BDS_steer(Group, PSind, beam0, model, Handles, Gind, 0.02);
  IPdat = IP_meas( Handles, PSind, beam0, model, [1 1] );
  s_disp(ind)=IPdat.Disp_y;
  s_coup(ind)=IPdat.sigma(2,3);
end
PS(PSind.Quad(Group.SQuad_id(5))).SetPt=p_init; stat=PSTrim(PSind.Quad(Group.SQuad_id(5)));
[q,dq,chisq,rms]=noplot_polyfit(pmoves,s_disp,1,1); ddispy_pm=q(2);
[q,dq,chisq,rms]=noplot_polyfit(pmoves,s_coup,1,1); dcoup_pm=q(2);
M=[ddispy_sm ddispy_pm;
   dcoup_sm  dcoup_pm];
y_knob{1}=M\[1;0]; y_knob{2}=M\[0;1];
y_knob{1}=y_knob{1}./abs(y_knob{1}(1)); y_knob{2}=y_knob{2}./abs(y_knob{2}(2));

% y_knob = -k2(s_list(end))*BEAMLINE{Group.Sext.dB.ClusterList(5).index(1)}.L*2*eta_s(end)*...
%     sqrt(model.Twiss.betay(Group.Sext.dB.ClusterList(5).index(1))*model.Twiss.betay(model.ip_ind))*...
%     sin((model.Twiss.nuy(model.ip_ind)-model.Twiss.nuy(Group.Sext.dB.ClusterList(5).index(1)))*2*pi)*2;

VERB=v_init;

return