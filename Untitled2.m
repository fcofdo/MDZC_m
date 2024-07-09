
clear
clc
Xm=500;
OpDis='rural';
kst=2;
xkm=Xm.*0.001;
[SYout]=SIGY(Xm,OpDis,kst);
[SZout]=SIGZ(Xm,OpDis,kst);
hs=30;
uref=4.3;
zref=10;
p=.15;
[us]=WSadj(hs,uref,zref,p)