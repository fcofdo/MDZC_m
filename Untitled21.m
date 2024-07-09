

angulo =0:11.25:360;
angulo(angulo>360) = angulo(angulo>360) - 360;
angulo(angulo<0)   = angulo(angulo(angulo<0)) + 360;
angulo=angulo';
tr=angulo/57.2958;
WDSIN=sin(tr);
WDCOS=cos(tr);
N=29;
XSourc=5000;
YSourc=5000;
a=0:100:10000;
[YRec,XRec] = meshgrid(a);
XRec=10000-XRec;

% DOWNWIND DISTANCE
X=(XRec-XSourc)*WDCOS(N) +(YRec-YSourc)*WDSIN(N);
% X=-((XRec-XSourc)*WDSIN(N)+(YRec-YSourc)*WDCOS(N)) ;
% Y=(YSourc-XRec)*WDCOS(N)+(XSourc-YRec)*WDSIN(N) ;
%  CROSSWIND DISTANCE
 Y=-((XRec-XSourc)*WDSIN(N)+(YRec-YSourc)*WDCOS(N));
% X=(YSourc-YRec)*WDSIN(N)+(XSourc-XRec)*WDCOS(N);
 X(X<0) = 0;
  Y(Y<0) = 0;
    
%  X = -((XR-XS)*WDSIN + (YR-YS)*WDCOS)
%       Y =   (XR-XS)*WDCOS - (YR-YS)*WDSIN
