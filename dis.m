
angulo =meteorologia(:,5);
% angulo =0:11.25:360;
angulo(angulo>360) = angulo(angulo>360) - 360;
angulo(angulo<0)   = angulo(angulo(angulo<0)) + 360;
tr=angulo/57.2958;
WDSIN=sin(tr);
WDCOS=cos(tr);
N=1;
XSourc=5000;
YSourc=5000;
a=0:100:10000;
[YRec,XRec] = meshgrid(a);
XRec=10000-XRec;
for i=1:length(angulo(:,1))
    % DOWNWIND DISTANCE
    X(:,:,i)=(XRec-XSourc)*WDCOS(i) +(YRec-YSourc)*WDSIN(i);
    % X=-((XRec-XSourc)*WDSIN(N)+(YRec-YSourc)*WDCOS(N)) ;
    % Y=(YSourc-XRec)*WDCOS(N)+(XSourc-YRec)*WDSIN(N) ;
    %  CROSSWIND DISTANCE
    Y(:,:,i)=-((XRec-XSourc)*WDSIN(i)+(YRec-YSourc)*WDCOS(i));
    % X=(YSourc-YRec)*WDSIN(N)+(XSourc-XRec)*WDCOS(N);
    X(:,:,i)=(X(:,:,i)<0) = 0;
    Y(:,:,i)=(Y(:,:,i)<0) = 0;
end
%  X = -((XR-XS)*WDSIN + (YR-YS)*WDCOS)
%       Y =   (XR-XS)*WDCOS - (YR-YS)*WDSIN
