function [tt,K,N,L,P,E,twss,orbt,S]=xtfft2mat5(fname)
%
% [tt,K,N,L,P,E,twss,orbt,S]=xtfft2mat5(fname);
%
% Outputs:
%
%   tt   = run title
%   K    = element keyword
%   N    = element name
%   L    = element length
%   P    = element parameter
%   E    = energy
%   twss = twiss (mux,betx,alfx,dx,dpx,muy,bety,alfy,dy,dpy)
%   orbt = orbit (x,px,y,py,t,pt)
%   S    = suml

% open the XTFF file

  fid=fopen(fname);
  if (fid==-1)
    error(['  Failed to open ',fname])
  end

% read in the header ... check that XTFF file is a twiss file

  line=fgetl(fid);
  xtff=line(9:16);
  if (~strcmp(xtff,'   TWISS'))
    error(['  Unexpected XTFF type (',xtff,') encountered ... abort'])
  end

% read in the run title

  tt=deblank(fgetl(fid));

% read in the INITIAL data

  line=fgetl(fid);
  K=line(1:4);
  N=line(5:20);
  L=str2double(line(21:32));
  p1=str2double(line(33:48));
  p2=str2double(line(49:64));
  p3=str2double(line(65:80));
  E=str2double(line(113:130));
  line=fgetl(fid);
  p4=str2double(line(1:16));
  p5=str2double(line(17:32));
  p6=str2double(line(33:48));
  p7=str2double(line(49:64));
  p8=str2double(line(65:80));
  line=fgetl(fid);
  t1=str2double(line(1:16));
  t2=str2double(line(17:32));
  t3=str2double(line(33:48));
  t4=str2double(line(49:64));
  t5=str2double(line(65:80));
  line=fgetl(fid);
  t6=str2double(line(1:16));
  t7=str2double(line(17:32));
  t8=str2double(line(33:48));
  t9=str2double(line(49:64));
  t10=str2double(line(65:80));
  line=fgetl(fid);
  orb1=str2double(line(1:16));
  orb2=str2double(line(17:32));
  orb3=str2double(line(33:48));
  orb4=str2double(line(49:64));
  S=str2double(line(65:80));
  P=[p1,p2,p3,p4,p5,p6,p7,p8];
  twss=[t3,t2,t1,t4,t5,t8,t7,t6,t9,t10];
  orbt=[orb1,orb2,orb3,orb4];

% read in the data ... break at end of the file

  while(1)
    line=fgetl(fid);
    if (length(line)==48),break,end
    K=[K;line(1:4)];
    N=[N;line(5:20)];
    L=[L;str2double(line(21:32))];
    p1=str2double(line(33:48));
    p2=str2double(line(49:64));
    p3=str2double(line(65:80));
    E=[E;str2double(line(113:130))];
    line=fgetl(fid);
    p4=str2double(line(1:16));
    p5=str2double(line(17:32));
    p6=str2double(line(33:48));
    p7=str2double(line(49:64));
    p8=str2double(line(65:80));
    line=fgetl(fid);
    t1=str2double(line(1:16));
    t2=str2double(line(17:32));
    t3=str2double(line(33:48));
    t4=str2double(line(49:64));
    t5=str2double(line(65:80));
    line=fgetl(fid);
    t6=str2double(line(1:16));
    t7=str2double(line(17:32));
    t8=str2double(line(33:48));
    t9=str2double(line(49:64));
    t10=str2double(line(65:80));
    line=fgetl(fid);
    orb1=str2double(line(1:16));
    orb2=str2double(line(17:32));
    orb3=str2double(line(33:48));
    orb4=str2double(line(49:64));
    S=[S;str2double(line(65:80))];
    P=[P;p1,p2,p3,p4,p5,p6,p7,p8];
    twss=[twss;t3,t2,t1,t4,t5,t8,t7,t6,t9,t10];
    orbt=[orbt;orb1,orb2,orb3,orb4];
  end

% close the XTFF file

  fclose(fid);
  
end