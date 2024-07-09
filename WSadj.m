function [us]=WSadj(hs,uref,zref,p)
%***********************************************************************
%                 WSADJ Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Adjusts Wind Speed from Anemometer Height to Stack Height

%        PROGRAMMER: Roger Brode, Jeff Wang

%        INPUTS:  Arrays of Source Parameters
%                 Meteorological Variables for One Hour
%                 Wind Speed Profile Exponents (Default or User-defined)

%        OUTPUTS: Stack Top Wind Speed, US

%        CALLED FROM:   PCALC
%                       VCALC
%                       ACALC
%***********************************************************************

%     Adjust Wind Speed -- Assume Wind Speed Constant Below 10 meters
if hs >= 10.0 &&  hs < 200.0
    us=uref.*(hs./zref).^p;
else if  hs >=  200.0
        us=uref.*(200./zref).^p;
    else if zref > 10.0 ;
            us = uref .*(10.0./zref).^p;
        else;
            us = uref;
        end
    end;
end %     Do Not Allow Stack Height Wind Speed < 1.0 m/s
us (us <1)=1;

 %subroutine wsadj

