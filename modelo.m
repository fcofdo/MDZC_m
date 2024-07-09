

clc;
clear;
% Colocar las coordenadas
XSourc=2500;
YSourc=5000;
Q=45;
hs=30;
zref=12;
OpDis='rural';
Vs=6.5;
D=0.80;
Ts=150;
Xmax=10000;
Ymax=10000;
N=25;
%cargar el archivo meteorologico
[Uadj,kst,angulo,hl,Ta]=meteo(hs,zref,OpDis);

% encuentra las distancias
[X,Y]=Dist_XY(angulo,XSourc,YSourc,Xmax,N);

%introduce la estabilidad y la zona

[SZout1,SYout1,H]=Sigmas(X,Ta,Ts,D,Vs,Uadj,hs,kst,OpDis);

z=zeros(size(X(:,:,1)));

[Conc]=ConcSourPoint(X,Q,Uadj,Y,SZout1,SYout1,z,H,kst,hl);
ConProme(Conc,Xmax,Ymax,N,YSourc, XSourc);
grabamovie(X,Conc,Xmax,Ymax,N,YSourc, XSourc);



