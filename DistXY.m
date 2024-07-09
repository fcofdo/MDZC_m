

function [XYdist]=DistXY(XR,XSRC,YR,YSRC,angulo)

%     Calculate Downwind (X) and Crosswind (Y) Distances
angulo(angulo>360) = angulo(angulo>360) - 360;
angulo(angulo<0)   = angulo(angulo(angulo<0)) + 360;
tr=angulo/57.2958;
WDSIN=sin(tr);
 WDCOS=cos(tr);


    XR=0:100:10000; XR= XR';     
      YR=0:100:10000; YR= YR';  
      XSRC=384;
      YSRC=175;
         SPA=[];
%                function [SPA]=ARDIST(XR,XSRC,YR,YSRC,angulo)

% for i=1:length(sp(:,1))
    for n=1:length(XR(:,1))
    R(n,1)=(XR(n,1)-XSRC(1));
    S(n,1)=YR(n,1)-YSRC(1);
    end
% end
for i=1:length(sp(:,1))
    for n=1:length(recp(:,1))
        for mi=1:length(angulo(:,1))
    X(i,n,mi)=S(i,n)*cost(mi,1)+R(i,n)*sint(mi,1);
    Y(i,n,mi)=S(i,n)*sint(mi,1)+R(i,n)*cost(mi,1);
        end
    end
end% end

for I = 1:length(XR)
         XVERT=XR(I)
         YVERT=YR(I)
         SPA(I,1) = -( XVERT-XSRC)*WDSIN(1) + (YVERT-YSRC)*WDCOS(1);
         SPA(I,2) =   ( XVERT-XSRC)*WDCOS(1) - (YVERT-YSRC)*WDSIN(1);
%          SPA(I,3)=DSQRT(SPA(I,1)*SPA(I,1) + SPA(I,2)*SPA(I,2))
end

