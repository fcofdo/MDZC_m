

function [Uadj,kst,angulo,hl,Ta]=meteo(hs,zref,OpDis)
%
load 'meteorologia.dat'

puz1=[0.07;0.07;0.10;0.15;0.35;0.55] ;
puz2=[0.15;0.15;0.20;0.25;0.30;0.30] ;
% for nm=1:lengeth(meteorologia(:,1))
angulo =meteorologia(:,5);
% tr=angulo/57.2958;
% sint=sin(tr);
% cost=cos(tr);
Uref=meteorologia(:,6);
kst=meteorologia(:,8);
% DTHDZ=0;
hl=meteorologia(:,9);
Ta=meteorologia(:,7);
pwr=zeros(length(kst(:,1)),1);
kst ((kst >6)) =6;
kst ((kst <1)) =1;
if OpDis=='rural'
    for I=1:length(kst(:,1))
        if kst(I,1)==1;pwr(I,1)=puz1(1,1); end
        if kst(I,1)==2;pwr(I,1)=puz1(2,1); end
        if kst(I,1)==3;pwr(I,1)=puz1(3,1); end
        if kst(I,1)==4;pwr(I,1)=puz1(4,1); end
        if kst(I,1)==5;pwr(I,1)=puz1(5,1); end
        if kst(I,1)==6;pwr(I,1)=puz1(6,1); end
    end
else
    for I1=1:length(kst(:,1))
        if kst(I1,1)==1;pwr(I1,1)=puz2(1,1); end
        if kst(I1,1)==2;pwr(I1,1)=puz2(2,1); end
        if kst(I1,1)==3;pwr(I1,1)=puz2(3,1); end
        if kst(I1,1)==4;pwr(I1,1)=puz2(4,1); end
        if kst(I1,1)==5;pwr(I1,1)=puz2(5,1); end
        if kst(I1,1)==6;pwr(I1,1)=puz2(6,1); end
    end
end

[Uadj]=WSadj(hs,Uref,zref,pwr);

angulo ((angulo >360)) = angulo ((angulo >360)) - 360;
angulo ((angulo <0))   = angulo ((angulo <0)) + 360;
% strColorMap = 'jet_20Intervals.txt';

% figure(1)
% wind_rose_ga(angulo,Uadj)
% %  set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperUnits', 'centimeters', 'PaperPosition', [0.25 0.25 28 19]);
% % % strOut = [strOutFileName1,'.bmp']
% %  strOut = ['rosa','.tif'];
% print('-dtiff','-r200','rosa')
% print( '-dtiff', '-r800', '-opengl', strOut);
% ;


