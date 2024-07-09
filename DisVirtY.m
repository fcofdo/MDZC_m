

function [XVY]=DisVirtY(SYO,KST)
% CALCULA LA DISTANCIA VIRTUAL NECESARIA PARA CONTABILIZAR LA DISPERSION
% INICIAL TRNSVERSAL
if KST==1
    XVY=(SYO/209.14)^0.890;
    if else KST==2
        XVY=(SYO/154.46)^0.902;
        if else KST==3
        XVY=(SYO/103.26)^0.917;
         if else KST==4
        XVY=(SYO/68.26)^0.919;
        if else KST==5
        XVY=(SYO/51.06)^0.921;
         if else KST==6
        XVY=(SYO/33.92)^0.919;
         end
        end
         end
        end
    end
end
return
