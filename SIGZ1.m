

function [SZout1]=SIGZ1(Xm,OpDis,kst)

%***********************************************************************
%      Calcula el valor de SIGMAZ - coeficiente de dispersion vertical
%        Programador: Francisco Garcia
%        Fecha:    Enero 2009
%        INPUTS:  Distancia eje X
%                 Clase de estabilidad
%                 Opcion de dispersión OpDis(Rural =1 o Urbana=2)
%       OUTPUTS: Vertical Dispersion Coefficient, SZout1

%        CALLED FROM:   PDIS
%                       VDIS
%                       ADIS
%                       SZENH
%                       DHPSS
%***********************************************************************
kst
km=Xm./1000

%     Determine Sigma-z Based on RURAL/URBAN, Stability Class, and Distance.
%     Stability Classes are Checked in the Order 4, 5, 6, 1, 2, 3
%     For Optimization, Since Neutral and Stable are Generally the Most
%     Frequent Classes.

if OpDis=='rural' %#ok<*STCMP>
    %        Retrieve Coefficients, A and B   ---   CALL SZCOEF
    if kst == 1
        if km <= 0.10
            A=122.8 ; b= 0.94470;
            SZout1=A*km^b;
        else if km <= 0.15
                A=158.080; b= 1.05420;
                SZout1=A*km^b;
            else if km <= 0.20
                    A=170.22; b= 1.09320;
                    SZout1=A*km^b;
                else if km <= 0.25
                        A=179.52;  b= 1.12620;
                        SZout1=A*km^b;
                    else if km <= 0.30
                            A=217.41;    b= 1.2644;
                            SZout1=A*km^b;
                        else if km <= 0.40
                                A=258.89; b= 1.4094;
                                SZout1=A*km^b;
                            else if km <= 0.50
                                    A=346.75 ;b= 1.72830;
                                    SZout1=A*km^b;
                                else
                                    A=453.85;b= 2.11660;
                                    SZout1=A*km^b;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if kst == 2
        if km <= 0.20
            A=90.673;b= 0.93198;
            SZout1=A*km^b;
        else if km <= 0.40
                A=98.483; b= 0.98332;
                SZout1=A*km^b;
            else
                A=109.3; b= 1.0971;
                SZout1=A*km^b;
            end
        end
    end
    
    if kst == 3
        A=61.141; b= 0.91465;
        SZout1=A*km^b;
    end
    if kst == 4
        if km <= .30
            A=34.459; b= 0.86974;
            SZout1=A*km^b;
        else if km <= 1.0
                A=32.093;  b= 0.81066;
                SZout1=A*km^b;
            else if km <= 3.0
                    A=32.093;b= 0.64403;
                    SZout1=A*km^b;
                else if km <= 10.
                        A=33.504; b= 0.60486;
                        SZout1=A*km^b;
                    else if km <= 30.
                            A=36.650;    b= 0.56589;
                            SZout1=A*km^b;
                        else    A=44.053; b= 0.51179;
                            SZout1=A*km^b;
                        end
                    end
                end
            end
        end
    end
    
    if kst == 5
        if km <= .10
            A=24.26; b= 0.83660;
            SZout1=A*km^b;
        else if km <= .30
                A=23.331; b= 0.81956;
                SZout1=A*km^b;
            else if km <= 1.0
                    A=21.628; b= 0.75660;
                    SZout1=A*km^b;
                else if km <= 2.0
                        A=21.628; b= 0.63077;
                        SZout1=A*km^b;
                    else if km <= 4.0
                            A=22.534; b= 0.57154;
                            SZout1=A*km^b;
                        else if km <= 10.
                                A=24.703; b= 0.50527;
                                SZout1=A*km^b;
                            else if km <= 20.
                                    A=26.97; b= 0.46713;
                                    SZout1=A*km^b;
                                else if km <= 40.
                                        A=35.42; b= 0.37615;
                                        SZout1=A*km^b;
                                    else
                                        A=47.618; b= 0.29592;
                                        SZout1=A*km^b;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if kst == 6
        if km <= .20
            A=15.209;    b= 0.81558;
            SZout1=A*km^b;
        else if km <= .70
                A=14.457;b= 0.78407;
                SZout1=A*km^b;
            else if km <= 1.0
                    A=13.953; b= 0.68465;
                    SZout1=A*km^b;
                else if km <= 2.0
                        A=13.953; b= 0.63227;
                        SZout1=A*km^b;
                    else if km <= 3.0
                            A=14.823;b= 0.54503;
                            SZout1=A*km^b;
                        else if km <= 7.0
                                A=16.187;b= 0.46490;
                                SZout1=A*km^b;
                            else if km <= 15.
                                    A=17.836;b= 0.41507;
                                    SZout1=A*km^b;
                                else if km <= 30.
                                        A=22.651; b= 0.32681;
                                        SZout1=A*km^b;
                                    else if km <= 60.
                                            A=27.074; b= 0.27436;
                                            SZout1=A*km^b;
                                        else
                                            A=34.219;   b= 0.21716;
                                            SZout1=A*km^b;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
else if OpDis=='urban';
        if kst == 4
            SZout1 = 140.*km./sqrt(1.+0.3.*km);
        else if kst >= 5 ;
                SZout1 = 80.*km./sqrt(1.+1.5.*km);
            else if kst <= 2 ;
                    SZout1 = 240.*km.*sqrt(1.+km);
                else if kst == 3 ;
                        SZout1 = 200.*km;
                    end;
                end;
            end
        end
    end
end
