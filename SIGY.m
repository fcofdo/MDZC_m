
function [SYyout]=SIGY(Xm,OpDis,kst);

%***********************************************************************
%                 SIGY Module of ISC2 Model
%        Calcula el valor de SIGMAY - coeficiente de dispersion horizontal
%        Programador: Francisco Garcia 
%        Fecha:    Enero 2009
%        INPUTS:  Distancia eje X
%                 Clase de estabilidad
%                 Opcion de dispersión OpDis(Rural =1 o Urbana=2)
%        OUTPUTS: Coefciciente e dispersion lateral, SYOUT
%        CALLED FROM:   PDIS
%                       VDIS
%                       ADIS
%                       SYENH
%                       DHPSS
%***********************************************************************
% clear global; clear functions;
xkm=Xm.*0.001;


%     Determine Sigma-y Based on RURAL/URBAN, Stability Class, and Distance.
%     Stability Classes are Checked in the Order 4, 5, 6, 1, 2, 3
%     For Optimization, Since Neutral and Stable are Generally the Most
%     Frequent Classes.
dtorad=0.017453293;

if OpDis=='rural'
    if kst == 4
        th =(8.3330 - 0.72382.*log(xkm)) .* dtorad;
    else if(kst == 5) ;
            th =(6.25 - 0.54287.*log(xkm)) .* dtorad;
        else if(kst == 6) ;
                th =(4.1667 - 0.36191.*log(xkm)) .* dtorad;
            else if(kst == 1) ;
                    th =(24.1667 - 2.5334.*log(xkm)) .* dtorad;
                else if(kst == 2) ;
                        th =(18.333 - 1.8096.*log(xkm)) .* dtorad;
                    else if(kst == 3) ;
                            th =(12.5 - 1.0857.*log(xkm)) .* dtorad;
                        end;
                    end
                end
            end
        end
    end
            
%        465.11628 = 1000. (m/km) / 2.15

SYyout = 465.11628 .* xkm .* tan(th);
if OpDis=='urban';
    if kst == 4
        SYyout = 160..*xkm./sqrt(1.+0.4.*xkm);
    else if kst >= 5 ;
            SYyout = 110..*xkm./sqrt(1.+0.4.*xkm);
        else if kst <= 2 ;
                SYyout = 320..*xkm./sqrt(1.+0.4.*xkm);
            else if kst == 3 ;
                    SYyout = 220..*xkm./sqrt(1.+0.4.*xkm);
                end;
            end
        end;
    end
end
end