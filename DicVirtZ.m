

function [XVZ]=DicVirtZ(SZinic,kst)
[A,B]=CoefXVZ(SZinic,kst);
XVZ=(SZinic/A)^B;

function [a,b]=CoefXVZ(SZinic,kst)
if kst == 1
    if SZinic<= 12.95
        a = 122.8 ; b = 1.0585;
    else if SZinic <= 21.40
            a = 158.080; b = 0.9486;
        else if SZinic <= 29.3
                a = 170.22; b = 0.9147;
            else if SZinic <= 37.67
                    a = 179.52;  b = 0.8879;
                else if SZinic <= 47.44
                        a = 217.41;    b = 0.7909;
                    else if SZinic <= 71.16
                            a = 258.89; b = 0.7095;
                        else if SZinic <=104.65
                                a = 346.75 ;b = 0.5786;
                            else
                                a = 453.85;b =0.4725;
                            end
                        end
                    end
                end
            end
        end
    end
end

if kst == 2
    if SZinic <= 20.23
        a = 90.673;b =1.073;
    else if SZinic <= 40
            a = 98.483; b = 1.017;
        else
            a = 109.3; b =0.9115;
        end
    end
end

if kst == 3
    a = 61.141; b = 1.0933;
end

if kst == 4
    if SZinic <= 12.09
        a = 34.459; b = 1.1498;
    else if SZinic <= 32.09
            a = 32.093;  b = 1.2336;
        else if SZinic <= 65.12
                a = 32.093;b = 1.5527;
            else if SZinic <=  134.9
                    a = 33.504; b =1.6533;
                else if SZinic <= 251.2
                        a = 36.650;    b = 1.7671;
                    else    a = 44.053; b = 1.9539;
                    end
                end
            end
        end
    end
end

if kst == 5
    if SZinic <= 3.534
        a = 24.26; b = 1.1953;
    else if SZinic <= 8.698
            a = 23.331; b = 1.2202;
        else if SZinic <= 21.628
                a = 21.628; b = 1.3217;
            else if SZinic <= 33.489
                    a = 21.628; b = 1.58547;
                else if SZinic <= 49.767
                        a = 22.534; b = 1.7497;
                    else if SZinic <= 79.07
                            a = 24.703; b = 1.9791;
                        else if SZinic <= 109.3
                                a = 26.97; b = 2.1407;
                            else if SZinic <= 141.86
                                    a = 35.42; b = 2.6585;
                                else
                                    a = 47.618; b = 31.3793;
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
    if SZinic <= 4.093
        a = 15.209;    b = 1.2261;
    else if SZinic <= 10.93
            a = 14.457;b = 1.2754;
        else if SZinic <= 13.953
                a = 13.953; b = 1.4606;
            else if SZinic <= 21.627
                    a = 13.953; b = 1.5816;
                else if SZinic <= 26.976
                        a = 14.823;b = 1.8348;
                    else if SZinic <= 40
                            a = 16.187;b = 2.151;
                        else if SZinic <= 54.89
                                a = 17.836;b = 2.4092;
                            else if SZinic <= 68.84
                                    a = 22.651; b = 3.0599;
                                else if SZinic <= 83.25
                                        a = 27.074; b = 3.16448;
                                    else
                                        a = 34.219;   b = 4.6049;
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


    
