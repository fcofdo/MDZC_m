

function [SYinic, SZinic]=SIGMAini(OpDis,kst,Ta,Ts,D,Vs,U,h)
%
Xm=0;
[DISTF,H,Delth]=Fplume(Ta,Ts,D,Vs,U,h,kst);
SYinic=Delth/3.5;
[SZout]=SIGZ(Xm,OpDis,kst);
SZinic=(SZout^2+(Delth/3.5)^2)^0.5;