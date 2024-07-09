

function [SZout]=SIGZ(Xm,OpDis,kst)

%***********************************************************************
%      Calcula el valor de SIGMAZ - coeficiente de dispersion vertical
%        Programador: Francisco Garcia
%        Fecha:    Enero 2009
%        INPUTS:  Distancia eje X
%                 Clase de estabilidad
%                 Opcion de dispersión OpDis(Rural =1 o Urbana=2)
%       OUTPUTS: Vertical Dispersion Coefficient, SZOUT

%        CALLED FROM:   PDIS
%                       VDIS
%                       ADIS
%                       SZENH
%                       DHPSS
%***********************************************************************

xkm =Xm/1000;

%     Determine Sigma-z Based on RURAL/URBAN, Stability Class, and Distance.
%     Stability Classes are Checked in the Order 4, 5, 6, 1, 2, 3
%     For Optimization, Since Neutral and Stable are Generally the Most
%     Frequent Classes.

if OpDis=='rural' %#ok<*STCMP>
    %        Retrieve Coefficients, A and B   ---   CALL SZCOEF
    [A,b]=SZcoef(xkm,kst);
    SZout=A*xkm^b;
else if OpDis=='urban';
        if kst == 4
            SZout = 140.*xkm./sqrt(1.+0.3.*xkm);
        else if kst >= 5 ;
                SZout = 80.*xkm./sqrt(1.+1.5.*xkm);
            else if kst <= 2 ;
                    SZout = 240.*xkm.*sqrt(1.+xkm);
                else if kst == 3 ;
                        SZout = 200.*xkm;
                    end;
                end;
            end
        end
    end
end
