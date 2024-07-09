
function [DistFx]=DistF(Ta,Ts,D,Vs,U,kst)
Ts=Ts+273.15;
for i=1:length(Ta)
    if Ta(i,1)<=0
        Ta(i,1)=293.15;
    end
end

F=2.4516375*Vs*(D^2)*(Ts-Ta)/Ta;

if kst<=4
    if F>=55
        
        XST=34.*F^(0.4);
    else
       
        XST=14.*F^(0.625);
    end
    DistFx=3.5*XST;
else if kst>=5
        if kst==5
            DTHDZ=0.02;
        else
            DTHDZ=0.035;
        end
        S=9.80616*DTHDZ/Ta;
        
        DistFx=2.0715*U/(S^(0.5));
    end
end
DistFx=DistFx/1000;


