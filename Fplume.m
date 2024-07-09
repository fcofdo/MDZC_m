
function [H]=Fplume(Ta,Ts,D,Vs,U,h,kst)
Ts=Ts+273.15;
for i=1:length(Ta)
    if Ta(i,1)<=0
        Ta(i,1)=293.15;
    end
end
if Vs<15*U
    hprim=h+2*((Vs/U)-1.5);
else hprim=h;
end
% Flujo volumetrico de gas (m3/s)
% Vf=0.785398*Vs*D*D;
F=2.4516375*Vs*(D^2)*(Ts-Ta)/Ta;

if kst<=4
    if F>=55
        Delth=(38.71*F^(3/5))/U;
        XST=34.*F^(0.4);
    else
        Delth=(21.425*F^(3/4))/U;
        XST=14.*F^(0.625);
    end
    DISTF=3.5*XST;
else if kst>=5
        if kst==5
            DTHDZ=0.02;
        else
            DTHDZ=0.035;
        end
        S=9.80616*DTHDZ/Ta;
        Delth=2.6*(F/(U*S))^(1/3);
        DISTF=2.0715*U/(S^(0.5));
    end
end
H=hprim+Delth;
DISTF=DISTF/1000;


