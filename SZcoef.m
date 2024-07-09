

function [A,b]=SZcoef(xkm,kst)

if kst == 1
    if xkm <= 0.10
        A=122.8 ; b= 0.94470;
    else if xkm <= 0.15
            A=158.080; b= 1.05420;
        else if xkm <= 0.20
                A=170.22; b= 1.09320;
            else if xkm <= 0.25
                    A=179.52;  b= 1.12620;
                else if xkm <= 0.30
                        A=217.41;    b= 1.2644;
                    else if xkm <= 0.40
                            A=258.89; b= 1.4094;
                        else if xkm <= 0.50
                                A=346.75 ;b= 1.72830;
                            else
                                A=453.85;b= 2.11660;
                            end
                        end
                    end
                end
            end
        end
    end
end

if kst == 2
    if xkm <= 0.20
        A=90.673;b= 0.93198;
    else if xkm <= 0.40
            A=98.483; b= 0.98332;
        else
            A=109.3; b= 1.0971;
        end
    end
end

if kst == 3
    A=61.141; b= 0.91465;
end
if kst == 4
    if xkm <= .30
        A=34.459; b= 0.86974;
    else if xkm <= 1.0
            A=32.093;  b= 0.81066;
        else if xkm <= 3.0
                A=32.093;b= 0.64403;
            else if xkm <= 10.
                    A=33.504; b= 0.60486;
                else if xkm <= 30.
                        A=36.650;    b= 0.56589;
                    else    A=44.053; b= 0.51179;
                    end
                end
            end
        end
    end
end

if kst == 5
    if xkm <= .10
        A=24.26; b= 0.83660;
    else if xkm <= .30
            A=23.331; b= 0.81956;
        else if xkm <= 1.0
                A=21.628; b= 0.75660;
            else if xkm <= 2.0
                    A=21.628; b= 0.63077;
                else if xkm <= 4.0
                        A=22.534; b= 0.57154;
                    else if xkm <= 10.
                            A=24.703; b= 0.50527;
                        else if xkm <= 20.
                                A=26.97; b= 0.46713;
                            else if xkm <= 40.
                                    A=35.42; b= 0.37615;
                                else
                                    A=47.618; b= 0.29592;
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
    if xkm <= .20
        A=15.209;    b= 0.81558;
    else if xkm <= .70
            A=14.457;b= 0.78407;
        else if xkm <= 1.0
                A=13.953; b= 0.68465;
            else if xkm <= 2.0
                    A=13.953; b= 0.63227;
                else if xkm <= 3.0
                        A=14.823;b= 0.54503;
                    else if xkm <= 7.0
                            A=16.187;b= 0.46490;
                        else if xkm <= 15.
                                A=17.836;b= 0.41507;
                            else if xkm <= 30.
                                    A=22.651; b= 0.32681;
                                else if xkm <= 60.
                                        A=27.074; b= 0.27436;
                                    else
                                        A=34.219;   b= 0.21716;
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




    
