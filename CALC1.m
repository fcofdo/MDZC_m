function calc(varargin)

% Date: 2012-07-28  Time: 13:43:37

%***********************************************************************
%                 CALC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Controls Flow and Processing of CALCulation Modules

%        PROGRAMMER: Roger Brode, Jeff Wang

%        DATE:    March 2, 1992

%        MODIFIED:   To move check for AREA source first for optimization.
%                    R. W. Brode, PES - 5/12/99

%        MODIFIED:   To add call for new source type of OPENPIT.
%                    R. W. Brode, PES - 9/30/94

%        INPUTS:  Arrays of Source Parameters
%                 Arrays of Receptor Locations
%                 Meteorological Variables for One Hour

%        OUTPUTS: Array of 1-hr CONC or DEPOS Values for Each Source/Receptor

%        CALLED FROM:   HRLOOP
%***********************************************************************

%     Variable Declarations
% use main1;
clear global; clear functions;

persistent modnam ; 

if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'CALC';
path = 'CN';

%     Begin Source LOOP
for isrc = 1: numsrc;
if(strcmp(deblank(srctyp(isrc)),deblank('area')))
%           Calculate Area Source Values for Rectangles  ---   CALL ACALC
acalc;
elseif(strcmp(deblank(srctyp(isrc)),deblank('point'))) ;
%           Calculate Point Source Values                ---   CALL PCALC
pcalc;
elseif(strcmp(deblank(srctyp(isrc)),deblank('volume'))) ;
%           Calculate Volume Source Values               ---   CALL VCALC
vcalc;
elseif(strcmp(deblank(srctyp(isrc)),deblank('areapoly'))) ;
%           Calculate Area Source Values for Polygons    ---   CALL ACALC
acalc;
elseif(strcmp(deblank(srctyp(isrc)),deblank('areacirc'))) ;
%           Calculate Area Source Values for Circles     ---   CALL ACALC
acalc;
elseif(strcmp(deblank(srctyp(isrc)),deblank('openpit'))) ;
%           Calculate OpenPit Source Values              ---   CALL OCALC
ocalc;
end;
end; isrc = numsrc+1;
%     end Source LOOP

return;
end %subroutine calc

function pcalc(varargin)
%***********************************************************************
%                 PCALC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates concentration or deposition values
%                 for POINT sources

%        PROGRAMMER: Roger Brode, Jeff Wang

%        MODIFIED BY D. Strimaitis, SRC (for WetDry DEPOSITION)
%        MODIFIED BY D. Strimaitis, SRC (for COMPLEX I -Intermediate
%                                        Terrain Processing)

%        DATE:    December 15, 1993

%        MODIFIED:   To skip receptor if it is more than 80km from source
%                    for TOXICS option.
%                    R. W. Brode, PES - 02/19/99

%        MODIFIED:   To allow use with EVENT processing.
%                    R. W. Brode, PES - 12/2/98

%        MODIFIED:   To add call for new source type of OPENPIT.
%                    R. W. Brode, PES - 9/30/94

%        MODIFIED BY D. Strimaitis, SRC (for Dry DEPOSITION)
%        (DATE:    February 15, 1993)

%        INPUTS:  Source Parameters for Specific Source
%                 Arrays of Receptor Locations
%                 Meteorological Variables for One Hour

%        OUTPUTS: 1-hr CONC or DEPOS Values for Each Receptor for
%                 Particular Source

%        CALLED FROM:   CALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent i modnam ; 
qs=[];ievent=[];irec=[];x=[];dhp=[];heflat=[];zelev=[];he=[];distr=[];dhpcmp=[];hecomp=[];hecmp1=[];corr=[];sy=[];sz=[];xy=[];xz=[];sbid=[];szcmp1=[];xzcmp1=[];sbcmp1=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(i), i=0; end;

%     Variable Initializations
modnam = 'PCALC';

%     Set the Source Variables for This Source              ---   CALL SETSRC
setsrc;

%     Apply Variable Emission Rate and Unit Factors         ---   CALL EMFACT
[qs]=emfact(qs);

%     Set Deposition Variables for this Source
if(ldpart || lwpart || ldgas)
%        Calculate Deposition Velocities for this Source    ---   CALL VDP
vdp;
end;

%     Calculate Scavenging Ratios for this Source           ---   CALL SCAVRAT
if(lwpart || lwgas)
scavrat;
end;

if((qtk ~= 0.0) &&(stable ||(hs <= zi) || depos || wdep))
%        Calculate Buoyancy and Momentum Fluxes
if(ts < ta)
ts = ta;
end;
fb =(0.25./ts).*(vs.*ds.*ds).*g.*(ts-ta);
fm =(0.25./ts).*(vs.*ds.*ds).*vs.*ta;
%        Adjust Wind Speed to Stack Height                  ---   CALL WSADJ
wsadj;
%        Calculate Distance to Final Rise                   ---   CALL DISTF
distf;
%        Set Wake and Building Type Switches                ---   CALL WAKFLG
wakflg;
%        Initialize FSTREC Logical Switch for First Receptor of Loop
fstrec = true;
%        Initialize FSTREC Logical Switch for First CMP1 Receptor of Loop
fstcmp = true;
if(ldpart || lwpart)
%           Calculate Min Sigma-z for Settled Plume @ Surface --- CALL SETSZMN
setszmn;
end;

%        Begin Receptor LOOP
for irec = 1: numrec;
%           Calculate Down and Crosswind Distances          ---   CALL XYDIST
if(evonly)
[ievent]=xydist(ievent);
else;
[irec]=xydist(irec);
end;
if(abs(y) > 1.191754.*x)
%              Receptor is at least 50 deg. off the plume centerline
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
elseif(distr < 0.99 || distr < 3..*zlb) ;
%              Receptor Too Close to Source for Calculation
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
elseif(toxics && distr > 80000.) ;
%              Receptor is beyond 80km from source.
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
else;

%RWB           Modifications to Integrate COMPLEX1 Algorithms.
%RWB           Get Plume Heights and Check for INTERMEDIATE TERRAIN Regime
if(nocmpl)
[x,dhp,heflat]=pheff(x,dhp,heflat);
[heflat,zelev,he]=sterad(heflat,zelev,he);
simple = true;
interm = false;
complx = false;
%                 Set HECOMP = HEFLAT for later check versus ZI
hecomp = heflat;
elseif(nosmpl) ;
%                 Note: Use radial distance, DISTR, for COMPLEX1 plume height.
[distr,dhpcmp,hecomp]=pheffc(distr,dhpcmp,hecomp);
[hecomp,zelev,hecmp1,corr]=cterad(hecomp,zelev,hecmp1,corr);
complx = true;
interm = false;
simple = false;
%                 Set HEFLAT = HECOMP for later check versus ZI
heflat = hecomp;
else;
[x,dhp,heflat]=pheff(x,dhp,heflat);
[heflat,zelev,he]=sterad(heflat,zelev,he);
if(elev)
[distr,dhpcmp,hecomp]=pheffc(distr,dhpcmp,hecomp);
[hecomp,zelev,hecmp1,corr]=cterad(hecomp,zelev,hecmp1,corr);
%                    Set the Simple/Intermediate/Complex Terrain Flags
itset;
else;
simple = true;
hecomp = heflat;
end;
end;

if(stable ||(heflat <= zi) ||(hecomp <= zi) ||depos  ||  wdep)

if(simple || interm)
%                    Determine Simple Terrain Sigmas        ---   CALL PDIS
[x,sy,sz,xy,xz,sbid]=pdis(x,sy,sz,xy,xz,sbid);
end;
if(complx || interm)
%                    Determine Complex Terrain Sigmas       ---   CALL PDISC
[distr,szcmp1,xzcmp1,sbcmp1]=pdisc(distr,szcmp1,xzcmp1,sbcmp1);
end;

%                 Determine Deposition Correction Factors for Gases
if(lwgas)
%                    Initialize wet source depletion factor to unity.
wqcorg = 1.;
wqcorgc = 1.;
if(wdplete)
%                       Determine source depletion factor
%                       from wet removal (GASES)
if(simple || interm)
%                          Simple Terrain Model
wqcorg = exp(-gscvrt.*x./us);
end;
if(complx || interm)
%                          Complex Terrain Model - use radial distance
wqcorgc = exp(-gscvrt.*distr./us);
end;
end;
end;

%                 Apply Intermediate Terrain Logic
if(simple)
%                    Simple Terrain Model Only              ---   CALL PSIMPL
psimpl;
elseif(complx) ;
%                    Complex Terrain Model Only             ---   CALL PCOMPL
pcompl;
elseif(interm) ;
%                    Initialize simple and complex terrain holding variables
simcon = 0.0;
comcon = 0.0;
for ityp = 1: numtyp;
simpl(ityp) = 0.;
compl(ityp) = 0.;
if(wetscim)
simpld(ityp) = 0.;
end;
if(wetscim)
compld(ityp) = 0.;
end;
end; ityp = numtyp+1;
%                    Determine Which Model Predicts the Larger Conc.
%                    Save Simple Terrain Conc.           ---   CALL PSIMPL
psimpl;
for ityp = 1: numtyp;
simpl(ityp) = hrval(ityp);
if(wetscim)
simpld(ityp) = hrvald(ityp);
end;
end; ityp = numtyp+1;
%                    Save Complex Terrain Conc.          ---   CALL PCOMPL
pcompl;
for ityp = 1: numtyp;
compl(ityp) = hrval(ityp);
if(wetscim)
compld(ityp) = hrvald(ityp);
end;
end; ityp = numtyp+1;
%                    Report Result for Model that Produces the Larger
%                    Concentration
if(simcon >= comcon)
for ityp = 1: numtyp;
hrval(ityp) = simpl(ityp);
if(wetscim)
hrvald(ityp) = simpld(ityp);
end;
end; ityp = numtyp+1;
else;
for ityp = 1: numtyp;
hrval(ityp) = compl(ityp);
if(wetscim)
hrvald(ityp) = compld(ityp);
end;
end; ityp = numtyp+1;
end;
end;

%                 Sum HRVAL to AVEVAL and ANNVAL Arrays  ---   CALL SUMVAL
if(evonly)
ev_sumval;
else;
sumval;
end;

else;
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
end;

%              Write DEBUG Information related to Terrain and Removal
if(debug)
writef(fid_iounit,['%0.15g \n']);
writef(fid_iounit,['%s %0.15g %0.15g \n'], 'HOUR, RECEPTOR : ',ihour,irec);
writef(fid_iounit,['%s %0.15g \n'], 'PCALC: HRVAL(final) = ',hrval);
if(ldpart || lwpart)
writef(fid_iounit,['%s \n'], 'PCALC: Particle Removal --------');
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'WQCOR  = ',wqcor(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'WQCORC = ',wqcorc(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'DQCOR  = ',dqcor(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'DQCORC = ',dqcorc(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'PCORZD = ',pcorzd(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'PCORZDC= ',pcorzdc(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'PCORZR = ',pcorzr(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'PCORZRC= ',pcorzrc(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'SZCOR  = ',szcor(i)); end;
for i=(1):(npd), writef(fid_iounit,['%s %0.15g \n'], 'SZCORC = ',szcorc(i)); end;
end;
writef(fid_iounit,['%s \n'], 'PCALC: Gas Removal -------------');
writef(fid_iounit,['%s %0.15g %0.15g \n'], 'WQCORG, WQCORGC = ',wqcorg,wqcorgc);
writef(fid_iounit,['%s \n'], 'PCALC: Concentration -----------');
writef(fid_iounit,['%s %0.15g %0.15g \n'], 'SIMPL, COMPL    = ',simpl,compl);
end;

end;
end; irec = numrec+1;
%        end Receptor LOOP
end;

return;
end %subroutine pcalc


function itset(varargin)
%***********************************************************************
%                 ITSET Module of the ISC Short Term Model - Version 2

%        PURPOSE:    To set intermediate terrain variables, based on
%                    complex terrain plume height and terrain height.

%        PROGRAMMER: Roger W. Brode, PES, Inc.

%        DATE:       September 30, 1994

%        INPUTS:     HECOMP = Complex Terrain Plume Height, without
%                             Terrain Adjustment Factors (through COMMON)

%        OUTPUTS:    SIMPLE, COMPLX, INTERM = Intermediate terrain
%                    logical control variables (through COMMON)

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent modnam ; 

if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'ITSET';
simple = false;
interm = false;
complx = false;

if(hs >= hter)
simple = true;
elseif(hecomp > hter) ;
interm = true;
else;
complx = true;
end;

%     Write Special DEBUG values for IT results
if(debug)
writef(fid_iounit,['%0.15g \n']);
writef(fid_iounit,['%s %0.15g \n'], 'ITSET --- IT RESULTS, HOUR :',ihour);
if(simple)
writef(fid_iounit,['%s \n'], 'ITFLAG = SIMPLE');
elseif(interm) ;
writef(fid_iounit,['%s \n'], 'ITFLAG = INTERM');
elseif(complx) ;
writef(fid_iounit,['%s \n'], 'ITFLAG = COMPLX');
end;
writef(fid_iounit,['%s %0.15g \n'], 'HECOMP         = ',hecomp);
writef(fid_iounit,['%s %0.15g %0.15g %0.15g \n'], 'HS, ZS, ZELEV  = ',hs,zs,zelev);
end;

return;
end %subroutine itset

function vcalc(varargin)
%***********************************************************************
%                 VCALC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates concentration or deposition values
%                 for VOLUME sources

%        PROGRAMMER: Roger Brode, Jeff Wang

%        MODIFIED BY D. Strimaitis, SRC (for WetDry DEPOSITION)

%        DATE:    December 15, 1993

%        MODIFIED:   To skip receptor if it is more than 80km from source
%                    for TOXICS option.
%                    R. W. Brode, PES - 02/19/99

%        MODIFIED:   To allow use with EVENT processing.
%                    R. W. Brode, PES - 12/2/98

%        MODIFIED BY D. Strimaitis, SRC (for Dry DEPOSITION)
%        (DATE:    February 15, 1993)

%        MODIFIED BY R. Brode, PES, to initialize SBID - 7/15/94

%        INPUTS:  Source Parameters for Specific Source
%                 Arrays of Receptor Locations
%                 Meteorological Variables for One Hour

%        OUTPUTS: 1-hr CONC or DEPOS Values for Each Receptor for
%                 Particular Source

%        CALLED FROM:   CALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent modnam ; 
qs=[];ievent=[];irec=[];zelev=[];heflat=[];he=[];x=[];sy=[];sz=[];xy=[];xz=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'VCALC';

%     Set the Source Variables for This Source              ---   CALL SETSRC
setsrc;

%     Apply Variable Emission Rate and Unit Factors         ---   CALL EMFACT
[qs]=emfact(qs);

%     Set Deposition Variables for this Source
if(ldpart || lwpart || ldgas)
%        Calculate Deposition Velocities for this Source    ---   CALL VDP
vdp;
end;

%     Calculate Scavenging Ratios for this Source           ---   CALL SCAVRAT
if(lwpart || lwgas)
scavrat;
end;

if((qtk ~= 0.0) &&(stable ||(hs <= zi) || depos || wdep))
%        Adjust Wind Speed to Release Height                ---   CALL WSADJ
wsadj;
%        Calculate Effective Radius
xrad = 2.15.*syinit;
if(ldpart || lwpart)
%           Calculate Min Sigma-z for Settled Plume @ Surface --- CALL SETSZMN
setszmn;
end;
%        Initialize SBID to 0.0 for call to DEPCOR
sbid = 0.0;

%        Begin Receptor LOOP
for irec = 1: numrec;
%           Calculate Down and Crosswind Distances          ---   CALL XYDIST
if(evonly)
[ievent]=xydist(ievent);
else;
[irec]=xydist(irec);
end;
if(abs(y) > 1.191754.*x)
%              Receptor is at least 50 deg. off the plume centerline
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
elseif(distr <(xrad+0.99)) ;
%              Receptor Too Close to Source for Calculation
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
elseif((x-xrad) < 0.0) ;
%              Receptor Upwind of Downwind Edge
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
elseif(toxics && distr > 80000.) ;
%              Receptor is beyond 80km from source.
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
else;
%              Determine Effective Plume Height             ---   CALL VHEFF
[zelev,heflat,he]=vheff(zelev,heflat,he);
%              Determine Dispersion Parameters              ---   CALL VDIS
[x,sy,sz,xy,xz]=vdis(x,sy,sz,xy,xz);
if(lwgas)
%                 Initialize wet source depletion factor to unity.
wqcorg = 1.;
if(wdplete)
%                    Determine source depletion factor
%                    from wet removal (GASES)
wqcorg=exp(-gscvrt.*x./us);
end;
end;
%              Calculate Conc. or Depos. for Virtual Point Source
%              Using a Simple Terrain Model                 ---   CALL PSIMPL
psimpl;

%              Sum HRVAL to AVEVAL and ANNVAL Arrays        ---   CALL SUMVAL
if(evonly)
ev_sumval;
else;
sumval;
end;

end;
end; irec = numrec+1;
%        end Receptor LOOP
end;

return;
end %subroutine vcalc


function acalc(varargin)
%***********************************************************************
%                 ACALC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates concentration or deposition values
%                 for AREA sources utilizing an integrated line source.

%        PROGRAMMER: Jeff Wang, Roger Brode

%        MODIFIED BY D. Strimaitis, SRC (for WetDry DEPOSITION)

%        DATE:    November 8, 1993

%        MODIFIED:   To modify distance at which model switches to point
%                    source approximation for area sources under TOXICS
%                    option.
%                    R. W. Brode, PES - 02/04/2002

%        MODIFIED:   To incorporate optimizations for TOXICS option.
%                    R. W. Brode, PES - 02/19/99

%        MODIFIED:   To allow use with EVENT processing and to call
%                    new ARDIST routine instead of XYDIST.
%                    R. W. Brode, PES - 12/2/98

%        MODIFIED by YICHENG ZHUANG, SRC to combine version 93188 with
%                 version 93046 - 9/28/93

%        MODIFIED:   To incorporate numerical integration algorithm
%                    for AREA source - 7/7/93

%        MODIFIED BY D. Strimaitis, SRC (for DEPOSITION) - 2/15/93

%        MODIFIED BY R. Brode, PES, to initialize XZ, XY, and SBID - 7/15/94

%*       MODIFIED BY J. Hardikar, PES, to make consistent with the new
%*                   OPENPIT Source Methodology - 7/20/94


%        INPUTS:  Source Parameters for Specific Source
%                 Arrays of Receptor Locations
%                 Meteorological Variables for One Hour

%        OUTPUTS: Array of 1-hr CONC or DEPOS Values for Each Source/Receptor

%        CALLED FROM:   CALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent i lterr modnam ; 
qs=[];ievent=[];xdep=[];width=[];length=[];xmaxr=[];irec=[];sy=[];xy=[];sz=[];xz=[];vdep=[];vgrav=[];zrdep=[];zflag=[];he=[];zi=[];us=[];xs=[];ys=[];xr=[];yr=[];rural=[];urban=[];kst=[];sbid=[];szmin=[];zelev=[];zs=[];debug=[];iounit=[];srctyp=[];isrc=[];ltgrid=[];kurdat=[];dqcor=[];pcorzr=[];pcorzd=[];szcor=[];toxics=[];vdepg=[];dqcorg=[];pcorzrg=[];pcorzdg=[];szcorg=[];x=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(i), i=0; end;
real :: :: xdep, width, length, xmaxr, qtksav, xpoint;
if isempty(lterr), lterr=false; end;

%     Variable Initializations
modnam = 'ACALC';
%     Set LTERR to falsemlv to signal simple terrain call to DEPCOR.
lterr = false;

%     Set the Source Variables for This Source              ---   CALL SETSRC
setsrc;

%     Apply Variable Emission Rate and Unit Factors         ---   CALL EMFACT
[qs]=emfact(qs);

%     Set Deposition Variables for this Source
if(ldpart || lwpart || ldgas)
%        Calculate Deposition Velocities for this Source    ---   CALL VDP
vdp;
end;

%     Calculate Scavenging Ratios for this Source           ---   CALL SCAVRAT
if(lwpart || lwgas)
scavrat;
end;

if((qtk ~= 0.0) &&(stable ||(hs <= zi) || depos || wdep))
%        Adjust Wind Speed to Release Height                ---   CALL WSADJ
wsadj;
if(ldpart || lwpart)
%           Calculate Min Sigma-z for Settled Plume @ Surface --- CALL SETSZMN
setszmn;
end;
%*       Initialize XY and XZ to 0.0 (XZ is used in
%*       call to DEPCOR from PLUMEF)
xy = 0.0;
xz = 0.0;

%        Initialize SBID to 0.0 (for call to DEPCOR from PLUMEF)
sbid = 0.0;

%        Begin Receptor LOOP
for irec = 1: numrec;
%           Calculate Down and Crosswind Distances          ---   CALL ARDIST
if(evonly)
[ievent,xdep,width,length,xmaxr]=ardist(ievent,xdep,width,length,xmaxr);
else;
[irec,xdep,width,length,xmaxr]=ardist(irec,xdep,width,length,xmaxr);
end;

%           Check to see if receptor is upwind of area source.
if(xmaxr < 1.0)
continue;
end;

%           Check to see if receptor is more than 80km from source.
if(toxics && distr > 80000.)
continue;
end;

%           Check to see if receptor is beyond edge of plume laterally.
if((abs(y)-0.5.*width) > 0.)
[xmaxr,sy,xy]=adisy(xmaxr,sy,xy);
if((abs(y)-0.5.*width) >= 4..*sy)
continue;
end;
end;

he = hs;
heflat = he;
if(stable ||(heflat <= zi) || depos || wdep)
%              Determine whether area integral or 'virtual' point will be used
%              Calculate distance for switch to point source approximation
xpoint = 1.5.*length + vp_fact.*width;
if(~toxics || strcmp(deblank(srctyp(isrc)),deblank('areapoly')) ||(toxics && x < xpoint))
if(ardplete && ldpart)
%                    Determine deposition correction factors for particles
syinit = 0.0;
[xdep,sz,xz]=adisz(xdep,sz,xz);
%                    Loop over particle sizes
for i = 1: npd;
%                       Initialize wetdry source depletion factors, profile
%                       correction factors, and sigma-z settling correction
%                       factors to unity.
wqcor(i)  = 1.0;
dqcor(i)  = 1.0;
pcorzr(i) = 1.0;
pcorzd(i) = 1.0;
szcor(i)  = 1.0;
%                       Determine factors for depletion
%                       from dry removal                    ---   CALL DEPCOR
[ vdep(i),vgrav(i),zrdep,zflag,xdep,xz,he,zi,us,xs,ys,xr,yr, rural,urban,kst,sz,sbid,szmin(i),zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat,dqcor(i),pcorzr(i),pcorzd(i),szcor(i),toxics]=depcor( vdep(i),vgrav(i),zrdep,zflag,xdep,xz,he,zi,us,xs,ys,xr,yr, rural,urban,kst,sz,sbid,szmin(i),zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat,dqcor(i),pcorzr(i),pcorzd(i),szcor(i),toxics);
end; i = fix(npd+1);
elseif(ardplete && ldgas) ;
%                    Determine deposition correction factors for particles
syinit = 0.0;
[xdep,sz,xz]=adisz(xdep,sz,xz);
%                    Initialize source depletion factors to unity.
dqcorg  = 1.0;
pcorzrg = 1.0;
pcorzdg = 1.0;
szcorg  = 1.0;
wqcorg  = 1.0;
%                    Determine factors for depletion - note that
%                    plume ht adjustment for terrain is signalled
%                    by a local logical - LTERR
%                    Simple Terrain Model                   ---   CALL DEPCOR
[ vdepg,dumvar2,zrdep,zflag, xdep,xz,he,zi,us,xs,ys,xr,yr,rural,urban,kst,sz,sbid,dumvar19,zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat, dqcorg,pcorzrg,pcorzdg,szcorg,toxics]=depcor( vdepg,0.0,zrdep,zflag, xdep,xz,he,zi,us,xs,ys,xr,yr,rural,urban,kst,sz,sbid, 2..*zrdep,zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat, dqcorg,pcorzrg,pcorzdg,szcorg,toxics);
end;
for ityp = 1: numtyp;
%                    Calculate Area Source Integral         ---   CALL AREAIN
areain;
end; ityp = numtyp+1;
else;
%                 use point source approximation
%                 Save emissions per unit area and calculate total emissions
qtksav = qtk;
if(strcmp(deblank(srctyp(isrc)),deblank('area')))
qtk = qtk .* xinit .* yinit;
elseif(strcmp(deblank(srctyp(isrc)),deblank('areacirc'))) ;
qtk = qtk .* pi .* radius(isrc) .* radius(isrc);
end;
syinit = 0.0;
[x,sy,sz,xy,xz]=vdis(x,sy,sz,xy,xz);
psimpl;
qtk = qtksav;
end;

else;
%              Plume Is Above Mixing Height, ZI, and No Wet Deposition
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
end;

%           Sum HRVAL to AVEVAL and ANNVAL Arrays           ---   CALL SUMVAL
if(evonly)
ev_sumval;
else;
sumval;
end;

end; irec = numrec+1;
%        end Receptor LOOP
end;

return;
end %subroutine acalc


function [indx,xdep,width,length,xmaxrec]=ardist(indx,xdep,width,length,xmaxrec);
%***********************************************************************
%                 ARDIST Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Sets Receptor Variables and Calculates Downwind (X)
%                 and Crosswind (Y) Distances, Crosswind Width (WIDTH),
%                 Distance used for AREADPLT Option (XDEP), Maximum
%                 Downwind Distance by Vertex (XMAXREC), and
%                 Radial Distance from Source to Receptor (DISTR)

%        PROGRAMMER: Roger Brode, Jeff Wang

%        DATE:    March 2, 1992

%        INPUTS:  Source Location
%                 Arrays of Receptor Locations
%                 SIN and COS of Wind Direction FROM Which Wind
%                 is Blowing, WDSIN and WDCOS

%        OUTPUTS: Values of X, Y, and DISTR (m) [in MAIN1]
%                 XDEP (m)
%                 WIDTH (m)
%                 LENGTH (m)
%                 XMAXREC (m)

%        CALLED FROM:   ACALC
%                       OCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent i modnam ; 

if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(i), i=0; end;
real :: :: xsrc, ysrc, xminrec, yminrec, ymaxrec,;

%     Variable Initializations
modnam = 'ARDIST';

%     Set Receptor Coordinates, Terrain Elevation and Flagpole Heights
xr = axr(indx);
yr = ayr(indx);
zelev = azelev(indx);
zflag = azflag(indx);

xminrec =  9999999.;
xmaxrec = -9999999.;
yminrec =  9999999.;
ymaxrec = -9999999.;

%     Calculate Downwind (X) and Crosswind (Y) Distances for Each Vertex
for i = 1: nvert+1;
xsrc = xvert(i);
ysrc = yvert(i);
spa(i,1) = -((xr-xsrc).*wdsin +(yr-ysrc).*wdcos);
spa(i,2) =(xr-xsrc).*wdcos -(yr-ysrc).*wdsin;
xminrec = min(xminrec, spa(i,1));
xmaxrec = max(xmaxrec, spa(i,1));
yminrec = min(yminrec, spa(i,2));
ymaxrec = max(ymaxrec, spa(i,2));
end; i = fix(nvert+1+1);

%     Calculate crosswind width, WIDTH, and alongwind length, LENGTH
width  = ymaxrec - yminrec;
length = xmaxrec - xminrec;

%     Determine downwind distance to use for AREADPLT option, XDEP
if(xminrec >= 0.0)
xdep = xminrec + 0.333333 .* length;
else;
xdep = 0.333333 .* xmaxrec;
end;

xdep = max( 1.0, xdep );

%     Calculate Downwind (X) and Crosswind (Y) Distances from Center of Source
x = -((xr-xcntr).*wdsin +(yr-ycntr).*wdcos);
y =(xr-xcntr).*wdcos -(yr-ycntr).*wdsin;

%     Calculate Radial Distance from Center of Source
distr = sqrt(x.*x + y.*y);

return;
end %subroutine ardist

function ocalc(varargin)
%***********************************************************************
%                 OCALC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates concentration or deposition values
%                 for OPENPIT sources

%        PROGRAMMER: Jayant Hardikar, Roger Brode
%        ADAPTED FROM:  SUBROUTINE ACALC

%        DATE:    July 19, 1994

%        MODIFIED:   To incorporate optimizations for TOXICS option.
%                    R. W. Brode, PES - 02/19/99

%        MODIFIED:   To allow use with EVENT processing.
%                    R. W. Brode, PES - 12/2/98

%        MODIFIED:   To skip calculations if QPTOT = 0.0, avoiding
%                    zero divide error in SUB. AMFRAC.
%                    R. W. Brode, PES Inc., - 4/14/95

%        INPUTS:  Source Parameters for Specific Source
%                 Arrays of Receptor Locations
%                 Meteorological Variables for One Hour

%        OUTPUTS: Array of 1-hr CONC or DEPOS Values for Each Source/Receptor

%        CALLED FROM:   CALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent lterr modnam ; 
icat=[];qptot=[];qeff=[];xr=[];yr=[];xvm=[];yvm=[];inout=[];ievent=[];xdep=[];width=[];length=[];xmaxr=[];irec=[];sz=[];xz=[];vdep=[];i=[];vgrav=[];zrdep=[];zflag=[];he=[];zi=[];us=[];xs=[];ys=[];rural=[];urban=[];kst=[];sbid=[];szmin=[];zelev=[];zs=[];debug=[];iounit=[];srctyp=[];isrc=[];ltgrid=[];kurdat=[];dqcor=[];pcorzr=[];pcorzd=[];szcor=[];toxics=[];vdepg=[];dqcorg=[];pcorzrg=[];pcorzdg=[];szcorg=[];x=[];sy=[];xy=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
integer :: :: i, ii, icat, inout;
real :: :: qptot, xvm(5), yvm(5), xdep, width, length, xmaxr, qtksav;
if isempty(lterr), lterr=false; end;

%     Variable Initializations
modnam = 'OCALC';
%     Set LTERR to falsemlv to signal simple terrain call to DEPCOR.
lterr = false;

%     Set the Source Variables for This Source              ---   CALL SETSRC
setsrc;

%*    Initialize the Total Adjusted Emission Rate from
%*    All Particles
qptot = 0.0;

%*    Loop over Particle Size Categories
for icat = 1:npd;
%*       Calculate the Escape Fraction for Each Category    ---   CALL ESCAPE
[icat]=escape(icat);

%*       Adjust the Emission Rate for Each Category         ---   CALL ADJEMI
[icat,qptot]=adjemi(icat,qptot);

%*    end Loop Over Particle Size Categories
end; icat =npd+1;

%*    Skip Calculations if QPTOT = 0.0
if(qptot == 0.0)
go to 999;
end;

%*    Adjust the Mass Fractions for All the Particle
%*    Size Categories                                       ---   CALL AMFRAC
[qptot]=amfrac(qptot);

%*    Determine the AlongWind Length of the OPENPIT Source  ---   CALL LWIND
lwind;

%*    Calculate the Relative Depth of the OPENPIT Source    ---   CALL PDEPTH
pdepth;

%*    Calculate the Fractional Size of the
%*    Effective Pit Area                                    ---   CALL PTFRAC
ptfrac;


%*    WRITE DEBUG INFORMATION
if(debug)
writef(fid_iounit,['%0.15g \n']);
writef(fid_iounit,['%0.15g \n']);
writef(fid_iounit,['%s \n'], 'DETAILED INFORMATION ON THE OPENPIT SOURCE:');
writef(fid_iounit,['%0.15g \n']);
writef(fid_iounit,['%0.15g \n']);
end;

%*    Determine the Coordinates of the Effective Pit Area
%*    in Wind Direction Coordinate System                   ---   CALL PITEFF
piteff;

%*    Calculate the Emission Rate for the Effective
%*    Pit Area                                              ---   CALL PITEMI
[qptot]=pitemi(qptot);

%*    WRITE DEBUG INFORMATION
if(debug)
writef(fid_iounit,['%s \n'], 'OPENPIT PARTICLE CHARACTERISTICS:');
writef(fid_iounit,['%s \n'], '---------------------------------');
writef(fid_iounit,['%0.15g \n']);
for ii =( 1):( npd), writef(fid_iounit,[repmat(' ',1,1),'ESCAPE FRACTIONS= ',repmat(['%8.3f',repmat(' ',1,2)] ,1,10) ' \n'],efrac(ii)); end;
%format (1X,'ESCAPE FRACTIONS= ',10(f8.3,2X));
for ii =( 1):( npd), writef(fid_iounit,[repmat(' ',1,1),'ADJUSTED EMISSION RATES= ',repmat(['%8.3f',repmat(' ',1,2)] ,1,10) ' \n'],qpart(ii)); end;
%format (1X,'ADJUSTED EMISSION RATES= ',10(f8.3,2X));
for ii =( 1):( npd), writef(fid_iounit,[repmat(' ',1,1),'ADJUSTED MASS FRACTIONS= ',repmat(['%8.3f',repmat(' ',1,2)] ,1,10) ' \n'],phi(ii)); end;
%format (1X,'ADJUSTED MASS FRACTIONS= ',10(f8.3,2X));
writef(fid_iounit,['%s %0.15g \n'], 'EMISSION RATE OF EFFECTIVE PIT= ',qeff);
writef(fid_iounit,['%0.15g \n']);
end;


%     Set Deposition Variables for this Source
if(ldpart || lwpart || ldgas)
%        Calculate Deposition Velocities for this Source    ---   CALL VDP
vdp;
end;

%     Calculate Scavenging Ratios for this Source           ---   CALL SCAVRAT
if(lwpart || lwgas)
scavrat;
end;

%     Apply Variable Emission Rate and Unit Factors         ---   CALL EMFACT
[qeff]=emfact(qeff);

if((qtk ~= 0.0) &&(stable ||(hs <= zi) || depos || wdep))
%        Adjust Wind Speed to Release Height                ---   CALL WSADJ
wsadj;
if(ldpart || lwpart)
%           Calculate Min Sigma-z for Settled Plume @ Surface --- CALL SETSZMN
setszmn;
end;
%        Initialize XY and XZ to 0.0 (XZ is used in call to DEPCOR from PLUMEF)
xy = 0.0;
xz = 0.0;
%        Initialize SBID to 0.0 (for call to DEPCOR from PLUMEF)
sbid = 0.0;
%        Begin Receptor LOOP
for irec = 1: numrec;
%           Check for receptor located inside boundary of open pit source
for i = 1: nvert+1;
xvm(i) = axvert(i,isrc);
yvm(i) = ayvert(i,isrc);
end; i = nvert+1+1;
xr = axr(irec);
yr = ayr(irec);
[xr,yr,xvm,yvm,dumvar5,inout]=pnpoly(xr,yr,xvm,yvm,5,inout);
if(inout > 0)
%              Receptor is within boundary - skip to next receptor
continue;
end;

%           Calculate Down and Crosswind Distances          ---   CALL ARDIST
if(evonly)
[ievent,xdep,width,length,xmaxr]=ardist(ievent,xdep,width,length,xmaxr);
else;
[irec,xdep,width,length,xmaxr]=ardist(irec,xdep,width,length,xmaxr);
end;

%           Check it see if receptor is upwind of area source
if(xmaxr < 1.0)
continue;
end;

%           Check to see if receptor is more than 80km from source.
if(toxics && distr > 80000.)
continue;
end;

he = hs;
heflat = he;
if(stable ||(heflat <= zi) || depos || wdep)
%              Determine whether area integral or 'virtual' point will be used
if(~toxics ||(toxics && x < vp_fact.*width))
if(ardplete && ldpart)
%                    Determine deposition correction factors for particles
syinit = 0.0;
[xdep,sz,xz]=adisz(xdep,sz,xz);
%                    Loop over particle sizes
for i = 1: npd;
%                       Initialize wetdry source depletion factors, profile
%                       correction factors, and sigma-z settling correction
%                       factors to unity.
wqcor(i)  = 1.0;
dqcor(i)  = 1.0;
pcorzr(i) = 1.0;
pcorzd(i) = 1.0;
szcor(i)  = 1.0;
%                       Determine factors for depletion
%                       from dry removal                    ---   CALL DEPCOR
[ vdep(i),vgrav(i),zrdep,zflag,xdep,xz,he,zi,us,xs,ys,xr,yr, rural,urban,kst,sz,sbid,szmin(i),zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat,dqcor(i),pcorzr(i),pcorzd(i),szcor(i),toxics]=depcor( vdep(i),vgrav(i),zrdep,zflag,xdep,xz,he,zi,us,xs,ys,xr,yr, rural,urban,kst,sz,sbid,szmin(i),zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat,dqcor(i),pcorzr(i),pcorzd(i),szcor(i),toxics);
end; i = npd+1;
elseif(ardplete && ldgas) ;
%                    Determine deposition correction factors for particles
syinit = 0.0;
[xdep,sz,xz]=adisz(xdep,sz,xz);
%                    Initialize source depletion factors to unity.
dqcorg  = 1.0;
pcorzrg = 1.0;
pcorzdg = 1.0;
szcorg  = 1.0;
wqcorg  = 1.0;
%                    Determine factors for depletion - note that
%                    plume ht adjustment for terrain is signalled
%                    by a local logical - LTERR
%                    Simple Terrain Model                   ---   CALL DEPCOR
[ vdepg,dumvar2,zrdep,zflag, xdep,xz,he,zi,us,xs,ys,xr,yr,rural,urban,kst,sz,sbid,dumvar19,zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat, dqcorg,pcorzrg,pcorzdg,szcorg,toxics]=depcor( vdepg,0.0,zrdep,zflag, xdep,xz,he,zi,us,xs,ys,xr,yr,rural,urban,kst,sz,sbid, 2..*zrdep,zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat, dqcorg,pcorzrg,pcorzdg,szcorg,toxics);
end;
for ityp = 1: numtyp;
%                    Calculate Area Source Integral            ---   CALL AREAIN
areain;
end; ityp = numtyp+1;
else;
%                 use point source approximation
%                 Save emissions per unit area and calculate total emissions
qtksav = qtk;
qtk = qtk .* xinit .* yinit;
syinit = 0.0;
[x,sy,sz,xy,xz]=vdis(x,sy,sz,xy,xz);
psimpl;
qtk = qtksav;
end;
else;
%              Plume Is Above Mixing Height, ZI, and No Wet Deposition
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
end;

%           Sum HRVAL to AVEVAL and ANNVAL Arrays        ---   CALL SUMVAL
if(evonly)
ev_sumval;
else;
sumval;
end;

end; irec = numrec+1;
%        end Receptor LOOP
end;

return;
end %subroutine ocalc


function setsrc(varargin)
%***********************************************************************
%                 SETSRC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Sets the Source Parameters for a Particular Source

%        PROGRAMMER: Roger Brode, Jeff Wang

%        MODIFIED BY: D. Strimaitis, SRC (for WetDry DEPOSITION)

%        DATE:  November 8,1993

%        MODIFIED by Yicheng Zhuang, SRC to combine version 93188 with
%                 version 93046 - 9/28/93

%        MODIFIED:   To incorporate inputs for numerical integration
%                    algorithm for AREA source - 7/7/93

%        MODIFIED BY D. Strimaitis, SRC (for DEPOSITION) - 2/15/93

%        MODIFIED BY Jayant Hardikar,PES (for handling OPENPIT
%                    Source - 7/19/94 , also modified AREA Source
%                    for Consistency with OPENPIT Source)

%        MODIFIED BY Roger Brode, PES (modified data structure for
%                    AXVERT and AYVERT for consistency with other
%                    2-D source arrays) - 8/15/95

%        INPUTS:  Source Parameters Arrays
%                 Source Index

%        OUTPUTS: Source Parameters for a Particular Source

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent j modnam ; 

if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(j), j=0; end;

%     Variable Initializations
modnam = 'SETSRC';

%     Assign The Values From Array Elements To Variables
if(strcmp(deblank(srctyp(isrc)),deblank('area')))
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
hs = ahs(isrc);

xinit = axinit(isrc);
yinit = ayinit(isrc);
angle = aangle(isrc);

szinit = aszini(isrc);
nvert = 4;

%        Store Vertices in Temporary Arrays
for ivert = 1: nvert+1;
xvert(ivert) = axvert(ivert,isrc);
yvert(ivert) = ayvert(ivert,isrc);
end; ivert = nvert+1+1;

xcntr = axcntr(isrc);
ycntr = aycntr(isrc);

elseif(strcmp(deblank(srctyp(isrc)),deblank('point'))) ;
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
hs = ahs(isrc);

ds = ads(isrc);
vs = avs(isrc);
ts = ats(isrc);

%        Check for Negative Stack Temperature, Used to Indicate Constant TS-TA
if(ts < 0.0)
ts = ta + abs(ts);
end;

if(ifvsec <= nsec)
dsbh = adsbh(ifvsec,isrc);
dsbw = adsbw(ifvsec,isrc);
if(idswak(ifvsec,isrc) == 0)
waklow = false;
elseif(idswak(ifvsec,isrc) == 1) ;
waklow = true;
end;
end;

elseif(strcmp(deblank(srctyp(isrc)),deblank('volume'))) ;
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
hs = ahs(isrc);

syinit = asyini(isrc);
szinit = aszini(isrc);

elseif(strcmp(deblank(srctyp(isrc)),deblank('area'))) ;
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
hs = ahs(isrc);

xinit = axinit(isrc);
yinit = ayinit(isrc);
angle = aangle(isrc);

szinit = aszini(isrc);
nvert = 4;

%        Store Vertices in Temporary Arrays
for ivert = 1: nvert+1;
xvert(ivert) = axvert(ivert,isrc);
yvert(ivert) = ayvert(ivert,isrc);
end; ivert = nvert+1+1;

xcntr = axcntr(isrc);
ycntr = aycntr(isrc);

elseif(strcmp(deblank(srctyp(isrc)),deblank('areapoly'))) ;
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
hs = ahs(isrc);

szinit = aszini(isrc);
nvert  = nverts(isrc);

%        Store Vertices in Temporary Arrays
for ivert = 1: nvert+1;
xvert(ivert) = axvert(ivert,isrc);
yvert(ivert) = ayvert(ivert,isrc);
end; ivert = nvert+1+1;

xcntr = axcntr(isrc);
ycntr = aycntr(isrc);

elseif(strcmp(deblank(srctyp(isrc)),deblank('areacirc'))) ;
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
hs = ahs(isrc);

szinit = aszini(isrc);
nvert  = nverts(isrc);

%        Store Vertices in Temporary Arrays
for ivert = 1: nvert+1;
xvert(ivert) = axvert(ivert,isrc);
yvert(ivert) = ayvert(ivert,isrc);
end; ivert = nvert+1+1;

xcntr = axcntr(isrc);
ycntr = aycntr(isrc);

elseif(strcmp(deblank(srctyp(isrc)),deblank('openpit'))) ;
xs = axs(isrc);
ys = ays(isrc);
zs = azs(isrc);
qs = aqs(isrc);
%        Set Emission Height of Effective Area, HS = 0.0
hs = 0.0;
%        Set Height of Emissions Above Base of Pit, EMIHGT
emihgt = ahs(isrc);
nvert = 4;

xinit = axinit(isrc);
yinit = ayinit(isrc);
angle = aangle(isrc);
palpha = aalpha(isrc);
pdeff  = apdeff(isrc);
szinit = aszini(isrc);
pitlen = max(xinit,yinit);
pitwid = min(xinit,yinit);

%        Store Vertices in Temporary Arrays
for ivert = 1: nvert+1;
xvert(ivert) = axvert(ivert,isrc);
yvert(ivert) = ayvert(ivert,isrc);
end; ivert = nvert+1+1;

xcntr = axcntr(isrc);
ycntr = aycntr(isrc);

end;

npd = inpd(isrc);
if(npd > 0)
for j = 1: npd;
pdiam(j) = apdiam(j,isrc);
phi(j) = aphi(j,isrc);
pdens(j) = apdens(j,isrc);
vgrav(j) = avgrav(j,isrc);
tstop(j) = atstop(j,isrc);
sc(j) = asc(j,isrc);
pscav(j,1) = apsliq(j,isrc);
pscav(j,2) = apsice(j,isrc);
end; j = fix(npd+1);
end;

%     Transfer Gas Wet Scavenging Coeff. (1:liquid, 2:frozen)
gscav(1) = agscav(1,isrc);
gscav(2) = agscav(2,isrc);

return;
end %subroutine setsrc

function [xarg,dhpout,heout]=pheff(xarg,dhpout,heout);
%***********************************************************************
%                 PHEFF Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Effective Plume Height for POINT Sources (m)

%        PROGRAMMER: Roger Brode, Jeff Wang

%        DATE:    March 2, 1992

%        MODIFIED:   To remove terrain adjustment to separate subroutine,
%                    and to use calling arguments
%                    R.W. Brode, PES, Inc. - 9/30/94

%        INPUTS:  Arrays of Source Parameters
%                 Logical Wake Flags
%                 Meteorological Variables for One Hour
%                 Wind Speed Adjusted to Stack Height
%                 Downwind Distance
%                 Terrain Elevation of Receptor

%        OUTPUTS: Plume Height (HEOUT) without Terrain Adjustment

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 
dhf=[];xf=[];
no type, intent(out)                     :: heout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
real :: :: heout, hsprim;

%     Variable Initializations
modnam = 'PHEFF';

%     Calculate Plume Height Without Terrain Adjustment
if((~ wake) &&(~ grdris))
%        Calculate Final Rise for First Receptor Only
if(fstrec)
fstrec = false;
hsp = hsprim(us,vs,hs,ds);
%           Calculate Final Rise, DHF                       ---   CALL DELH
[dhf]=delh(dhf);
end;
if(nostd)
heout = hs + dhf;
else;
heout = hsp + dhf;
end;
if(~ nobid)
%           Calculate Gradual Plume Rise for Use in BID Calculation
if(xarg < xf)
%              Calculate Gradual Rise, DHPOUT               ---   CALL DHPHS
[xarg,dhf,dhpout]=dhphs(xarg,dhf,dhpout);
else;
dhpout = dhf;
end;
else;
dhpout = dhf;
end;
elseif(wake && wakess) ;
%        Calculate Final Rise for First Receptor Only
if(fstrec)
fstrec = false;
%           Calculate Final Rise (at X=XF), DHF             ---   CALL DHPSS
[xf,dhpout]=dhpss(xf,dhpout);
dhf    = dhpout;
end;
if(xarg < xf)
%           Calculate Gradual Rise, DHP                     ---   CALL DHPSS
[xarg,dhpout]=dhpss(xarg,dhpout);
else;
dhpout = dhf;
end;
heout = hs + dhpout;
else;
%RWB       if ((WAKE .AND. (.NOT. WAKESS)) .OR.
%RWB          ((.NOT. WAKE) .AND. GRDRIS)) then
%        Calculate Final Rise for First Receptor Only
if(fstrec)
fstrec = false;
hsp = hsprim(us,vs,hs,ds);
%           Calculate Final Rise, DHF                       ---   CALL DELH
[dhf]=delh(dhf);
end;
if(xarg < xf)
%           Calculate Gradual Rise, DHP                     ---   CALL DHPHS
[xarg,dhf,dhpout]=dhphs(xarg,dhf,dhpout);
else;
dhpout = dhf;
end;
if(nostd)
heout = hs + dhpout;
else;
heout = hsp + dhpout;
end;
end;

return;
end %subroutine pheff

function [xarg,dhpout,heout]=pheffc(xarg,dhpout,heout);
%***********************************************************************
%                 PHEFFC Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Effective Plume Height for POINT Sources (m)
%                 in Complex Terrain

%        PROGRAMMER: Roger Brode

%        DATE:    September 30, 1994

%        INPUTS:  Arrays of Source Parameters
%                 Logical Wake Flags
%                 Meteorological Variables for One Hour
%                 Wind Speed Adjusted to Stack Height
%                 Downwind Distance
%                 Terrain Elevation of Receptor

%        OUTPUTS: Plume Height (HEOUT) without Terrain Adjustment

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 
dhfcmp=[];
no type, intent(out)                     :: heout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
real :: :: heout, hsprim;

%     Variable Initializations
modnam = 'PHEFFC';

if(fstcmp)
fstcmp = false;
%        This is the First Call for PHEFFC - Calculate HSP and DHFCMP
hsp = hsprim(us,vs,hs,ds);
[dhfcmp]=delh(dhfcmp);
end;

if(xarg < xf)
%        Distance is less than distance to final rise - Calculate gradual rise
[xarg,dhfcmp,dhpout]=dhphs(xarg,dhfcmp,dhpout);
else;
%        Set gradual rise = final rise for XARG > XF
dhpout = dhfcmp;
end;

%     Check for stack-tip downwash option
if(nostd)
heout = hs + dhpout;
else;
heout = hsp + dhpout;
end;

return;
end %subroutine pheffc

function [hearg,zarg,heout]=sterad(hearg,zarg,heout);
%***********************************************************************
%                 STERAD Module of the ISC Short Term Model - Version 2

%        PURPOSE: Adjusts Effective Plume Height for Simple Terrain Effects

%        PROGRAMMER: Roger Brode

%        DATE:    September 30, 1994

%        INPUTS:  HEARG = Flat terrain plume height
%                 ZARG  = Elevation of terrain

%        OUTPUTS: HEOUT = Effective plume height with terrain adjustment

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 

no type, intent(out)                     :: heout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
real :: :: heout, hterchop;

%     Variable Initializations
modnam = 'STERAD';

%     Adjust Plume Height for Elevated Terrain, Save Flat Terrain Value (HEFLAT)
%     For Later Comparison With Mixing Height
if(flat)
heout  = hearg;
elseif(elev) ;
%        Calculate Terrain Hgt Above Plant Grade (Chopped-off at Release Height)
hterchop = min( hs,(zarg - zs));
heout = hearg - hterchop;
end;

%     Don't Allow Effective Plume Height to be < 0.0
heout = max( 0.0, heout);

return;
end %subroutine sterad

function [hearg,zarg,heout,cout]=cterad(hearg,zarg,heout,cout);
%***********************************************************************
%                 CTERAD Module of the ISC Short Term Model - Version 2

%        PURPOSE: Adjusts Effective Plume Height for Complex Terrain Effects

%        PROGRAMMER: Roger Brode

%        DATE:    September 30, 1994

%        INPUTS:  HEARG = Flat terrain plume height
%                 ZARG  = Elevation of terrain

%        OUTPUTS: HEOUT = Effective plume height with terrain adjustment
%                 COUT  = Attenuation correction factor

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 

no type, intent(out)                     :: cout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'CTERAD';

%     Calculate Terrain Hgt Above Plant Grade
hter = zarg - zs;

%     Calculate COMPLEX1 Plume Height
heout = max((hearg.*tcf(kst)),(hearg-(1.0-tcf(kst)).*hter) );
heout = max( heout, zmin );

%     Calculate the Attentuation Correction Factor, COUT
if((unstab||neutrl) ||(hearg >=(hter+zflag)) )
cout = 1.0;
elseif((hter+zflag-hearg) >= 400.) ;
cout = 0.0;
else;
cout =(400. -(hter+zflag-hearg))./400.;
end;

return;
end %subroutine cterad

function [zarg,hefout,heout]=vheff(zarg,hefout,heout);
%***********************************************************************
%                 VHEFF Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Effective Plume Height for VOLUME Sources (m)

%        PROGRAMMER: Roger Brode, Jeff Wang

%        DATE:    March 2, 1992

%        INPUTS:  Arrays of Source Parameters
%                 Logical Wake Flags
%                 Meteorological Variables for One Hour
%                 Wind Speed Adjusted to Stack Height
%                 Downwind Distance
%                 Terrain Elevation of Receptor

%        OUTPUTS: Effective Plume Height (HE)

%        CALLED FROM:   VCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 

no type, intent(out)                     :: heout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
real :: :: heout, hterchop;

%     Variable Initializations
modnam = 'VHEFF';

%     Calculate Terrain Height Above Plant Grade (Chopped-off at Release Height)
if(flat)
hterchop = 0.0;
elseif(elev) ;
hterchop = min( hs,(zarg - zs));
end;

%     Calculate Effective Plume Height (No Rise) Adjusted for Terrain Height
heout = hs - hterchop;

%     Save Plume Height for Flat Terrain for Later Comparison to Mixing Height
hefout = hs;

return;
end %subroutine vheff

function [xarg,syout,szout,xyout,xzout,sbout]=pdis(xarg,syout,szout,xyout,xzout,sbout);
%***********************************************************************
%                 PDIS Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Dispersion Parameters for POINT Sources

%        PROGRAMMER: Roger Brode, Jeff Wang
%        MODIFIED BY D. Strimaitis, SRC (initialize SBID to 0.0)

%        DATE:    February 15, 1993

%        INPUTS:  Arrays of Source Parameters
%                 Logical Wake Flags
%                 Wake Plume Height, HEMWAK
%                 Meteorological Variables for One Hour
%                 Downwind Distance

%        OUTPUTS: Lateral and Vertical Dispersion Coefficients, SY and SZ

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 
syarg=[];szarg=[];dhp=[];
no type, intent(out)                     :: sbout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
real :: :: sbout, syarg, szarg;

%     Variable Initializations
modnam = 'PDIS';

if(~ wake)
%        Calculate Sigma-y from Curves                   ---   CALL SIGY
[xarg,syarg]=sigy(xarg,syarg);
%        Calculate Sigma-z from Curves                   ---   CALL SIGZ
[xarg,szarg]=sigz(xarg,szarg);
if(~ nobid)
%           Apply BID                                    ---   CALL BID
[dhp,syarg,szarg,syout,szout,sbout]=bid(dhp,syarg,szarg,syout,szout,sbout);
else;
sbout = 0.0;
syout = syarg;
szout = szarg;
end;
xyout = 0.0;
xzout = 0.0;
elseif(wake) ;
if(hemwak > 1.2.*dsbh)
%           Calculate Sigma-y from Curves                ---   CALL SIGY
[xarg,syarg]=sigy(xarg,syarg);
xyout = 0.0;
else;
%           Calculate Building Enhanced Sigma-y          ---   CALL SYENH
[xarg,syarg,xyout]=syenh(xarg,syarg,xyout);
end;
%        Calculate Building Enhanced Sigma-z             ---   CALL SZENH
[xarg,szarg,xzout]=szenh(xarg,szarg,xzout);
if((~ nobid) &&(~ wakess))
%           Apply BID                                    ---   CALL BID
[dhp,syarg,szarg,syout,szout,sbout]=bid(dhp,syarg,szarg,syout,szout,sbout);
else;
sbout = 0.0;
syout = syarg;
szout = szarg;
end;
end;

if(szout > 5000. && npd == 0)
szout = 5000.;
end;

return;
end %subroutine pdis

function [xarg,szout,xzout,sbout]=pdisc(xarg,szout,xzout,sbout);
%***********************************************************************
%                 PDISC Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Dispersion Parameters for POINT Sources

%        PROGRAMMER: Roger Brode

%        DATE:    March 2, 1992

%        INPUTS:  Arrays of Source Parameters
%                 Logical Wake Flags
%                 Wake Plume Height, HEMWAK
%                 Meteorological Variables for One Hour
%                 Downwind Distance

%        OUTPUTS: Lateral and Vertical Dispersion Coefficients, SY and SZ

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 
szarg=[];dhpcmp=[];syarg=[];syout=[];
no type, intent(out)                     :: sbout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
real :: :: sbout, syout, syarg, szarg;

%     Variable Initializations
modnam = 'PDISC';

%     Calculate Sigma-z from Curves Using Radial Distance   ---   CALL SIGZ
[xarg,szarg]=sigz(xarg,szarg);
syarg = 0.0;

if(~ nobid)
%        Apply BID                                          ---   CALL BID
[dhpcmp,syarg,szarg,syout,szout,sbout]=bid(dhpcmp,syarg,szarg,syout,szout,sbout);
else;
sbout = 0.0;
szout = szarg;
syout = 0.0;
end;
xzout = 0.0;

if(szout > 5000. && npd == 0)
szout = 5000.;
end;

return;
end %subroutine pdisc

function [xarg,syout,szout,xyout,xzout]=vdis(xarg,syout,szout,xyout,xzout);
%***********************************************************************
%                 VDIS Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Dispersion Parameters for VOLUME Sources

%        PROGRAMMER: Roger Brode, Jeff Wang

%        DATE:    March 2, 1992

%        MODIFIED:   To remove input distance calling argument for XVZ
%                    R. W. Brode, PES, Inc. - 12/29/97

%        INPUTS:  Arrays of Source Parameters
%                 Meteorological Variables for One Hour
%                 Downwind Distance

%        OUTPUTS: Lateral and Vertical Dispersion Coefficients

%        CALLED FROM:   VCALC
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 

no type, intent(in out)                  :: xyout;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'VDIS';

%     Calculate Lateral Virtual Distance                 ---   CALL XVY
[xyout]=xvy(xyout);
%     Calculate Sigma-y from Curves for X+XY             ---   CALL SIGY
[dumvar1,syout]=sigy(xarg+xyout,syout);
%     Calculate Vertical Virtual Distance                ---   CALL XVZ
[xzout]=xvz(xzout);
%     Calculate Sigma-z from Curves for X+XZ             ---   CALL SIGZ
[dumvar1,szout]=sigz(xarg+xzout,szout);

if(szout > 5000. && npd == 0)
szout = 5000.;
end;

return;
end %subroutine vdis

function [xarg,syout,xyout]=adisy(xarg,syout,xyout);
%***********************************************************************
%                 ADISY Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Lateral Dispersion Parameters for AREA Sources

%        PROGRAMMER: Roger Brode, PES, Inc.

%        DATE:    December 14, 1998

%        INPUTS:  Downwind Distance (m), XARG
%                 Meteorological Variables for One Hour

%        OUTPUTS: Lateral Dispersion Coefficient (m), SYOUT
%                 Lateral Virtual Distance (m), XYOUT

%        CALLED FROM:   PLUMEF, PWIDTH
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 

no type, intent(in out)                  :: xarg;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'ADISY';

%     Calculate Sigma-y from Curves for XARG                ---   CALL SIGY
[xarg,syout]=sigy(xarg,syout);
syout = max(syout,0.0001);
xyout = 0.0;

return;
end %subroutine adisy

function [xarg,szout,xzout]=adisz(xarg,szout,xzout);
%***********************************************************************
%                 ADISZ Module of the ISC Short Term Model - Version 2

%        PURPOSE: Calculates Vertical Dispersion Parameters for AREA Sources

%        PROGRAMMER: Roger Brode, PES, Inc.

%        DATE:    December 14, 1998

%        INPUTS:  Downwind Distance (m), XARG
%                 Meteorological Variables for One Hour

%        OUTPUTS: Vertical Dispersion Coefficient (m), SZOUT
%                 Vertical Virtual Distance (m), XZOUT

%        CALLED FROM:   PLUMEF
%***********************************************************************

%     Variable Declarations
% use main1;

persistent modnam ; 

no type, intent(in out)                  :: xarg;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'ADISZ';

%     Calculate Sigma-z from Curves for XARG                ---   CALL SIGZ
[xarg,szout]=sigz(xarg,szout);
szout = max(szout,0.0001);
xzout = 0.0;

%     Add Initial Dispersion for OPENPIT Sources
if(szinit > 0.0)
szout = sqrt(szout.*szout + szinit.*szinit);
end;

if(szout > 5000. && npd == 0)
szout = 5000.;
end;

return;
end %subroutine adisz

function psimpl(varargin)
%***********************************************************************
%               PSIMPL Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Hourly Concentration or Deposition
%                 value for POINT Sources
%                 Using Gaussian Plume Equation for Simple Terrain

%                 (Replaces PCHI and PDEP)

%           NOTE: Particle settling is treated as a 'tilted plume'
%                 until the centerline reaches the surface.  Thereafter
%                 the centroid height of the plume continues to be
%                 modified by gravity.  This process is simulated by
%                 altering the sigma-z for each particle-size.  Hence,
%                 sigma-z is now a function of particle-size.

%        PROGRAMMER: D. Strimaitis, SRC

%        DATE:    December 15, 1993

%        MODIFIED:   To compare YTERM to -18.0 rather than EXPLIM,
%                    equivalent to 6*SY.  R.W. Brode, PES, Inc. - 02/19/99

%        MODIFIED:   To call PDEP for call to SUB. DEPCOR; to use
%                    modified SUB. VERT.  R.W. Brode, PES, Inc. - 9/30/94

%        INPUTS:

%        OUTPUTS: HRVAL, Concentration or Deposition for Particular
%                 Source/Receptor Combination

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent j modnam wdonly ; 
x=[];he=[];sz=[];a0=[];zflag=[];v=[];ityp=[];zrdep=[];vd=[];vsimp=[];hesetl=[];szadj=[];vj=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(j), j=0; end;
real :: :: yterm, vsimp, a0, szadj, hv, adj, vj, dryflux, wetflux,vterm, vsimpd, vtmp;
if isempty(wdonly), wdonly=false; end;

%     Variable Initializations
modnam = 'PSIMPL';
wdonly = false;

if((unstab || neutrl) && heflat > zi)
%        Plume Is Above Mixing Height, ZI
if(depos || wdep)
%           Set WDONLY flag for Wet Deposition Only
wdonly = true;
else;
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
return;
end;
end;

yterm = -0.5.*(y.*y)./(sy.*sy);

%     If receptor is more than 6.*SY off the centerline (YTERM < -18.),
%     then skip calculation, otherwise continue.
if(yterm > -18.0)

if(npd == 0)
%           Determine Deposition Correction Factors for Gases
if(ldgas || lwgas)
[x, wdonly]=pdepg(x, wdonly);
end;
for ityp = 1: numtyp;
v(ityp) = 0.;
if(wetscim)
vdry(ityp) = 0.;
end;
end; ityp = numtyp+1;
vsimp = 0.0;
if(wetscim)
vsimpd = 0.0;
end;
ityp = 0;
a0  = -0.5./(sz.*sz);
adj = dqcorg .* wqcorg;
if(conc)
%              Concentration
ityp = ityp + 1;
if(wdonly)
%                 Plume is above mixing height so set CONC = 0.0
v(ityp) = 0.0;
else;
%                 Calculate Concentration Form of V         ---   CALL VERT
[he,sz,a0,zflag,v(ityp)]=vert(he,sz,a0,zflag,v(ityp));
v(ityp) = adj.*pcorzrg.*v(ityp)./sz;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
end;
if(depos || ddep)
if(wdonly)
%                 Plume is above mixing height so set DDEP = 0.0
dryflux = 0.0;
else;
%                 For Dry Deposition Complete Vertical Term is Needed
%                 Calculated at ZRDEP                       ---   CALL VERT
[he,sz,a0,zrdep,vd]=vert(he,sz,a0,zrdep,vd);
dryflux = adj.*pcorzdg.*vdepg.*vd./sz;
end;
end;
if(depos || wdep)
%              Calculate Wet Flux
%              For Wet Flux, Vertical Term is Integral of EXP terms
%              Over All z, so VJ/SZ=SQRT(2PI)
wetflux = adj.*gscvrt.*srt2pi;
end;
if(depos)
%              WetDry fluxes of particles are summed
ityp = ityp + 1;
v(ityp) = dryflux + wetflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
if(ddep)
%              Dry flux of particles
ityp = ityp + 1;
v(ityp) = dryflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
if(wdep)
%              Wet flux of particles
ityp = ityp + 1;
v(ityp) = wetflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
if(~conc && interm)
if(wdonly)
%                 Plume is above mixing height so set CONC = 0.0
vsimp = 0.0;
else;
%                 For Concentration Complete Vertical Term is Needed for
%                 Each Particle Size Calculated at ZFLAG ---   CALL VERT
[he,sz,a0,zflag,vsimp]=vert(he,sz,a0,zflag,vsimp);
%                 Calculate Concentration for Intermediate Terrain Check
vsimp = adj.*pcorzrg.*vsimp./sz;
if(wetscim)
vsimpd = vsimp ./ wqcorg;
end;
end;
end;

else;
%           Determine Deposition Correction Factors for Particles
if(ldpart || lwpart)
[x, wdonly]=pdep(x, wdonly);
end;

%           Calculate the Vertical Term, V for particles
for ityp = 1: numtyp;
v(ityp) = 0.;
if(wetscim)
vdry(ityp) = 0.;
end;
end; ityp = numtyp+1;
vsimp = 0.;
if(wetscim)
vsimpd = 0.;
end;
for j = 1: npd;
ityp = 0;
%              Settling may alter SZ for the Jth particle plume
szadj = sz.*szcor(j);
a0 = -0.5./(szadj.*szadj);
%              Calculate Plume Tilt Due to Settling, HV
hv =(x./us) .* vgrav(j);
%              Calculate Settled Plume Height, HESETL
hesetl = he - hv;
%              Restrict settled height to be positive, so that the plume
%              does not settle below the surface -- this is the limit of
%              the tilted plume technique.
hesetl = max(0.0,hesetl);
%              Adjust Jth contribution by mass fraction and source
%              depletion
adj = phi(j) .* dqcor(j) .* wqcor(j);
if(conc)
%                 Concentration
ityp = ityp + 1;
if(wdonly)
%                    Plume is above mixing height so set CONC = 0.0
v(ityp) = 0.0;
else;
%                    For Concentration Complete Vertical Term is Needed for
%                    Each Particle Size Calculated at ZFLAG ---   CALL VERT
[hesetl,szadj,a0,zflag,vj]=vert(hesetl,szadj,a0,zflag,vj);
vtmp = adj.*pcorzr(j).*vj./szadj;
v(ityp) = v(ityp) + vtmp;
if(wetscim)
vdry(ityp) = vdry(ityp) +vtmp./wqcor(j);
end;
end;
end;
if(depos || ddep)
if(wdonly)
%                    Plume is above mixing height so set DDEP = 0.0
dryflux = 0.0;
else;
%                    For Dry Deposition Complete Vertical Term is Needed for
%                    Each Particle Size Calculated at ZRDEP ---   CALL VERT
[hesetl,szadj,a0,zrdep,vj]=vert(hesetl,szadj,a0,zrdep,vj);
%                    Calculate Dry Flux VJ/SZ
dryflux = adj.*pcorzd(j).*vdep(j).*vj./szadj;
end;
end;
if(depos || wdep)
%                 Calculate Wet Flux VJ/SZ --
%                 For Wet Flux, Vertical Term is Integral of EXP terms
%                 Over All z, so VJ/SZ=SQRT(2PI)
wetflux = adj.*pscvrt(j).*srt2pi;
end;
if(depos)
%                 WetDry fluxes of particles are summed
ityp = ityp + 1;
vtmp = dryflux + wetflux;
v(ityp) = v(ityp) + vtmp;
if(wetscim)
vdry(ityp) = vdry(ityp) + vtmp./wqcor(j);
end;
end;
if(ddep)
%                 Dry flux of particles
ityp = ityp + 1;
v(ityp) = v(ityp) + dryflux;
if(wetscim)
vdry(ityp) = vdry(ityp) +dryflux./wqcor(j);
end;
end;
if(wdep)
%                 Wet flux of particles
ityp = ityp + 1;
v(ityp) = v(ityp) + wetflux;
if(wetscim)
vdry(ityp) = vdry(ityp) +wetflux./wqcor(j);
end;
end;
if(~conc && interm)
if(wdonly)
%                    Plume is above mixing height so set CONC = 0.0
vsimp = 0.0;
else;
%                    For Concentration Complete Vertical Term is Needed for
%                    Each Particle Size Calculated at ZFLAG ---   CALL VERT
[hesetl,szadj,a0,zflag,vj]=vert(hesetl,szadj,a0,zflag,vj);
%                    Calculate Concentration for Intermediate Terrain Check
vtmp  = adj.*pcorzr(j).*vj./szadj;
vsimp = vsimp + vtmp;
if(wetscim)
vsimpd = vsimpd + vtmp ./ wqcor(j);
end;
end;
end;
end; j = fix(npd+1);
end;

%        Calculate the Decay Term, D                     ---   CALL DECAY
[x]=decay(x);

for ityp = 1: numtyp;
%           Complete VTERM (SZ already in denomenator of V)
vterm =(d.*v(ityp))./(twopi.*us.*sy);

%           Check for Possible Underflow Condition
if(vterm > 0.0 &&(log(vterm)+yterm) > explim)
hrval(ityp) = qtk .* emifac(ityp) .* vterm .* exp(yterm);
else;
hrval(ityp) = 0.0;
end;
if(wetscim)
%              Repeat the above calculations for HRVALD array

vterm =(d.*vdry(ityp))./(twopi.*us.*sy);
if(vterm > 0.0 &&(log(vterm)+yterm) > explim)
hrvald(ityp) = qtk .* emifac(ityp) .* vterm .* exp(yterm);
else;
hrvald(ityp) = 0.0;
end;
end;
end; ityp = numtyp+1;

if(~conc && interm)
%           Calculate Concentration for Simple Terrain
%           Complete VTERM (SZ already in denomenator of V)
vterm =(d.*vsimp)./(twopi.*us.*sy);

%           Check for Possible Underflow Condition
if(vterm > 0.0 &&(log(vterm)+yterm) > explim)
simcon = qtk .* emicon .* vterm .* exp(yterm);
else;
simcon = 0.0;
end;
elseif(conc && interm) ;
simcon = hrval(1);
end;

else;
%        Lateral Term is 0.0; Set HRVAL's to 0.0
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
simcon = 0.0;
end;

if(debug)
%        Print Out Debugging Information                    ---   CALL DEBOUT
debout;
end;

return;
end %subroutine psimpl

function pcompl(varargin)
%***********************************************************************
%               PCOMPL Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Hourly Concentration or Deposition
%                 value for POINT Sources
%                 Using Gaussian Plume Equation for Complex Terrain

%                 (Replaces PCHI and PDEP)

%           NOTE: Particle settling is treated as a 'tilted plume'
%                 until the centerline reaches the surface.  Thereafter
%                 the centroid height of the plume continues to be
%                 modified by gravity.  This process is simulated by
%                 altering the sigma-z for each particle-size.  Hence,
%                 sigma-z is now a function of particle-size.

%        PROGRAMMER: D. Strimaitis, SRC

%        DATE:    December 15, 1993

%        MODIFIED:   To use fully integrated COMPLEX1 algorithms rather
%                    than calls to CMP1.  R.W. Brode, PES, Inc. - 9/30/94

%        INPUTS:

%        OUTPUTS: HRVAL, Concentration or Deposition for Particular
%                 Source/Receptor Combination

%        CALLED FROM:   PCALC
%***********************************************************************

%     Variable Declarations
% use main1;
persistent j modnam wdonly ; 
hecmp1=[];szcmp1=[];a0=[];zflag=[];v=[];ityp=[];zrdep=[];vd=[];vcomp=[];hesetl=[];zelev=[];corrj=[];szadj=[];vj=[];x=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(j), j=0; end;
real :: :: vcomp, a0, szadj, hv, corrj, adj, vj, dryflux, wetflux,wcmp1, vcompd, vtmp;
if isempty(wdonly), wdonly=false; end;

%     Variable Initializations
modnam = 'PCOMPL';
wdonly = false;

if((unstab || neutrl) && hecomp > zi)
%        Plume Is Above Mixing Height, ZI
if(depos || wdep)
%           Set WDONLY flag for Wet Deposition Only
wdonly = true;
else;
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
return;
end;
end;

if(abs(y) <= x.*0.19891 && corr > 0.0)
%        Receptor is inside of sector and Plume is < 400m Below Receptor

if(npd == 0)
%           Determine Deposition Correction Factors for Gases
if(ldgas || lwgas)
[wdonly]=pdepgc(wdonly);
end;
for ityp = 1: numtyp;
v(ityp) = 0.;
if(wetscim)
vdry(ityp) = 0.;
end;
end; ityp = numtyp+1;
vcomp = 0.0;
if(wetscim)
vcompd = 0.0;
end;
ityp = 0;
a0 = -0.5./(szcmp1.*szcmp1);
adj = dqcorgc .* wqcorgc;
if(conc)
%              Concentration
ityp = ityp + 1;
if(wdonly)
%                 Plume is above mixing height so set CONC = 0.0
v(ityp) = 0.0;
if(wetscim)
vdry(ityp) = 0.0;
end;
else;
%                 For Concentration, Complete Vertical Term is Needed for
%                 Each Particle Size Calculated at ZFLAG ---   CALL VERT
[hecmp1,szcmp1,a0,zflag,v(ityp)]=vert(hecmp1,szcmp1,a0,zflag,v(ityp));
v(ityp) = corr.*adj.*pcorzrgc.*v(ityp)./szcmp1;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorgc;
end;
end;
end;
if(depos || ddep)
if(wdonly)
%                 Plume is above mixing height so set DRYFLUX = 0.0
dryflux = 0.0;
else;
%                 For Dry Deposition, Complete Vertical Term is Needed for
%                 Each Particle Size Calculated at ZRDEP ---   CALL VERT
[hecmp1,szcmp1,a0,zrdep,vd]=vert(hecmp1,szcmp1,a0,zrdep,vd);
%                 Calculate Dry Flux VJ/SZ
dryflux = corr.*adj.*pcorzdgc.*vdepg.*vd./szcmp1;
end;
end;
if(depos || wdep)
%              Calculate Wet Flux VJ/SZ --
%              For Wet Flux, Vertical Term is Integral of EXP terms
%              Over All z, so VJ/SZ=SQRT(2PI)
wetflux = corr.*adj.*gscvrt.*srt2pi;
end;
if(depos)
%              WetDry fluxes of particles are summed
ityp = ityp + 1;
v(ityp) = dryflux + wetflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorgc;
end;
end;
if(ddep)
%              Dry flux of particles
ityp = ityp + 1;
v(ityp) = dryflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorgc;
end;
end;
if(wdep)
%              Wet flux of particles
ityp = ityp + 1;
v(ityp) = wetflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorgc;
end;
end;
if(~conc && interm)
if(wdonly)
%                 Plume is above mixing height so set CONC = 0.0
vcomp = 0.0;
else;
%                 For Concentration, Complete Vertical Term is Needed for
%                 Each Particle Size Calculated at ZFLAG ---   CALL VERT
[hecmp1,szcmp1,a0,zflag,vcomp]=vert(hecmp1,szcmp1,a0,zflag,vcomp);
%                 Calculate Concentration for Intermediate Terrain Check
vcomp = corr.*adj.*pcorzrgc.*vcomp./szcmp1;
if(wetscim)
vcompd = vcomp ./ wqcorgc;
end;
end;
end;

else;
%           Determine Deposition Correction Factors for Particles
if(ldpart || lwpart)
[wdonly]=pdepc(wdonly);
end;
%           Calculate the Vertical Term, V for particles
for ityp = 1: numtyp;
v(ityp) = 0.;
if(wetscim)
vdry(ityp) = 0.;
end;
end; ityp = numtyp+1;
vcomp = 0.;
if(wetscim)
vcompd = 0.0;
end;
for j = 1: npd;
ityp = 0;
%              Settling may alter SZ for the Jth particle plume
szadj = szcmp1.*szcorc(j);
a0 = -0.5./(szadj.*szadj);
%              Calculate Plume Tilt Due to Settling, HV
hv =(x./us) .* vgrav(j);
%              Calculate Settled Plume Height, HESETL
hesetl = hecomp - hv;
%              Restrict settled height to be positive, so that the plume
%              does not settle below the surface -- this is the limit of
%              the tilted plume technique.
hesetl = max(0.0,hesetl);
%              Calculate Adjusted Plume Height and Attenuation Factor
%              for This Particle Category
[hesetl,zelev,hecmp1,corrj]=cterad(hesetl,zelev,hecmp1,corrj);
%              Adjust Jth contribution by mass fraction and source
%              depletion
adj = phi(j) .* dqcorc(j) .* wqcorc(j);
if(conc)
%                 Concentration
ityp = ityp + 1;
if(wdonly)
%                    Plume is above mixing height so set CONC = 0.0
v(ityp) = 0.0;
else;
%                    For Concentration, Complete Vertical Term is Needed for
%                    Each Particle Size Calculated at ZFLAG ---   CALL VERT
[hecmp1,szadj,a0,zflag,vj]=vert(hecmp1,szadj,a0,zflag,vj);
vtmp = corrj.*adj.*pcorzrc(j).*vj./szadj;
v(ityp) = v(ityp) + vtmp;
if(wetscim)
vdry(ityp) = vdry(ityp)+vtmp./wqcorc(j);
end;
end;
end;
if(depos || ddep)
if(wdonly)
%                    Plume is above mixing height so set DRYFLUX = 0.0
dryflux = 0.0;
else;
%                    For Dry Deposition, Complete Vertical Term is Needed for
%                    Each Particle Size Calculated at ZRDEP ---   CALL VERT
[hecmp1,szadj,a0,zrdep,vj]=vert(hecmp1,szadj,a0,zrdep,vj);
%                    Calculate Dry Flux VJ/SZ
dryflux = corrj.*adj.*pcorzdc(j).*vdep(j).*vj./szadj;
end;
end;
if(depos || wdep)
%                 Calculate Wet Flux VJ/SZ --
%                 For Wet Flux, Vertical Term is Integral of EXP terms
%                 Over All z, so VJ/SZ=SQRT(2PI)
wetflux = corrj.*adj.*pscvrt(j).*srt2pi;
end;
if(depos)
%                 WetDry fluxes of particles are summed
ityp = ityp + 1;
vtmp = dryflux + wetflux;
v(ityp) = v(ityp) + vtmp;
if(wetscim)
vdry(ityp) = vdry(ityp) + vtmp./wqcorc(j);
end;
end;
if(ddep)
%                 Dry flux of particles
ityp = ityp + 1;
v(ityp) = v(ityp) + dryflux;
if(wetscim)
vdry(ityp) = vdry(ityp)+dryflux./wqcorc(j);
end;
end;
if(wdep)
%                 Wet flux of particles
ityp = ityp + 1;
v(ityp) = v(ityp) + wetflux;
if(wetscim)
vdry(ityp) = vdry(ityp)+wetflux./wqcorc(j);
end;
end;
if(~conc && interm)
if(wdonly)
%                    Plume is above mixing height so set CONC = 0.0
vcomp = 0.0;
else;
%                    For Concentration, Complete Vertical Term is Needed for
%                    Each Particle Size Calculated at ZFLAG ---   CALL VERT
[hecmp1,szadj,a0,zflag,vj]=vert(hecmp1,szadj,a0,zflag,vj);
%                    Calculate Concentration for Intermediate Terrain Check
vtmp  = corrj.*adj.*pcorzrc(j).*vj./szadj;
vcomp = vcomp + vtmp;
vcompd= vcompd+ vtmp./wqcorc(j);
end;
end;
end; j = fix(npd+1);
end;

%        Calculate the Decay Term, D                        ---   CALL DECAY
[x]=decay(x);

for ityp = 1: numtyp;
%           Calculate HRVAL for Sector Average in Complex Terrain
hrval(ityp) =(qtk.*emifac(ityp).*d.*v(ityp)) ./(srt2pi.*distr.*delthp.*us);
if(wetscim)
%              Repeat the above calculations for HRVALDHRVALJD arrays
hrvald(ityp) =(qtk.*emifac(ityp).*d.*vdry(ityp)) ./(srt2pi.*distr.*delthp.*us);
end;
end; ityp = numtyp+1;

if(~conc && interm)
%           Calculate Concentration for Sector Average in Complex Terrain
comcon =(qtk.*emicon.*d.*vcomp) ./(srt2pi.*distr.*delthp.*us);
elseif(conc && interm) ;
comcon = hrval(1);
end;

else;
%        Receptor is outside of sector or Plume is > 400m Below Receptor
for ityp = 1: numtyp;
hrval(ityp) = 0.0;
if(wetscim)
hrvald(ityp) = 0.0;
end;
end; ityp = numtyp+1;
comcon = 0.0;
end;

if(debug)
wcmp1 = delthp .* distr;
writef(fid_iounit,['%s \n'], 'PCOMPL ----------------------------------');
writef(fid_iounit,['%s %0.15g %0.15g \n'], 'Hour, Receptor     =',ihour,irec);
writef(fid_iounit,['%s \n'], '  ');
writef(fid_iounit,['%s %0.15g %0.15g \n'], 'QTK, D             =',qtk,d);
writef(fid_iounit,['%s %0.15g %0.15g %0.15g \n'], 'CORRJ, WCMP1, US   =',corrj,wcmp1,us);
writef(fid_iounit,['%s \n'], 'PCOMPL ----------------------------------');
end;

if(debug)
%        Print Out Debugging Information                    ---   CALL DEBOUT
%RWB         CALL DEBOUT
end;

return;
end %subroutine pcompl


function [xarg,rcz,rczd]=asimpl(xarg,rcz,rczd);
%***********************************************************************
%               ASIMPL Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Hourly Concentration or Deposition
%                 value for AREA Sources Using Numerical
%                 Integration Algorithm for Simple Terrain

%                 (Replaces ACHI and ADEP)

%           NOTE: Particle settling is treated as a 'tilted plume'
%                 until the centerline reaches the surface.  Thereafter
%                 the centroid height of the plume continues to be
%                 modified by gravity.  This process is simulated by
%                 altering the sigma-z for each particle-size.  Hence,
%                 sigma-z is now a function of particle-size.

%        PROGRAMMER: D. Strimaitis, SRC

%        DATE:    December 15, 1993

%        MODIFIED:   To correct problem if CONC DDEP and WDEP are selected
%                    for area sources (DDEP results were repeated for WDEP
%                    under previous version).  Also corrects test for
%                    WDONLY (wet deposition only) flag.
%                    R. W. Brode, PES, Inc. - 12/29/97

%        INPUTS:  Downwind Distance (m), XARG

%        OUTPUTS: Relative Vertical Component of Concentration or Deposition
%                 for A Unit Of Source/Receptor Combination, RCZ

%        CALLED FROM:   PLUMEF
%***********************************************************************

%     Variable Declarations
% use main1;

persistent j modnam sconc sddep sdepos swdep wdonly ; 
he=[];sz=[];a0=[];zflag=[];v=[];ityp=[];zrdep=[];vtmp=[];hesetl=[];szadj=[];vj=[];
no type, intent(out)                     :: rcz;
no type, intent(out)                     :: rczd;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(j), j=0; end;
real :: :: rcz, a0, szadj, hv, adj, vj, dryflux, wetflux;
real :: :: rczd, vtmp;
if isempty(sconc), sconc=false; end;
if isempty(sdepos), sdepos=false; end;
if isempty(sddep), sddep=false; end;
if isempty(swdep), swdep=false; end;
if isempty(wdonly), wdonly=false; end;

%     Variable Initializations
modnam = 'ASIMPL';
wdonly = false;

%     Determine appropriate output type for this ITYP, assign output type
%     logicals to local variables, and set others to false
if(strcmp(deblank(outtyp(ityp)),deblank('conc')))
sconc  = true;
sdepos = false;
sddep  = false;
swdep  = false;
elseif(strcmp(deblank(outtyp(ityp)),deblank('depos'))) ;
sconc  = false;
sdepos = true;
sddep  = false;
swdep  = false;
elseif(strcmp(deblank(outtyp(ityp)),deblank('ddep'))) ;
sconc  = false;
sdepos = false;
sddep  = true;
swdep  = false;
elseif(strcmp(deblank(outtyp(ityp)),deblank('wdep'))) ;
sconc  = false;
sdepos = false;
sddep  = false;
swdep  = true;
end;

if((unstab || neutrl) && heflat > zi)
%        Plume is above mixing height, ZI
if(sdepos || swdep)
%           Set WDONLY flag for Wet Deposition Only
wdonly = true;
else;
v(ityp) = 0.0;
rcz     = 0.0;
if(wetscim)
vdry(ityp) = 0.0;
end;
if(wetscim)
rczd = 0.0;
end;
return;
end;
end;

rcz  = 0.0;
rczd = 0.0;
if(xarg >= 1.0)
if(npd == 0)
%           Determine Deposition Correction Factors for Gases
if(~ardplete &&(ldgas || lwgas) )
[xarg, wdonly]=pdepg(xarg, wdonly);
end;
v(ityp) = 0.0;
vdry(ityp) = 0.0;
a0  = -0.5./(sz.*sz);
adj = dqcorg .* wqcorg;
%           Calculate the Vertical Term, V, for gases
if(sconc)
if(wdonly)
%                 Plume is above mixing height so set CONC = 0.0
v(ityp) = 0.0;
if(wetscim)
vdry(ityp) = 0.0;
end;
else;
%                 Calculate Concentration Form of V         ---   CALL VERT
[he,sz,a0,zflag,v(ityp)]=vert(he,sz,a0,zflag,v(ityp));
v(ityp) = adj.*pcorzrg.*v(ityp)./sz;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
end;
if(sdepos || sddep)
if(wdonly)
%                 Plume is above mixing height so set DDEP = 0.0
dryflux = 0.0;
else;
%                 For Dry Deposition Complete Vertical Term is Needed
%                 Calculated at ZRDEP                       ---   CALL VERT
[he,sz,a0,zrdep,vtmp]=vert(he,sz,a0,zrdep,vtmp);
dryflux = adj.*pcorzdg.*vdepg.*vtmp./sz;
end;
end;
if(sdepos || swdep)
%              Calculate Wet Flux
%              For Wet Flux, Vertical Term is Integral of EXP terms
%              Over All z, so VJ/SZ=SQRT(2PI)
wetflux = adj.*gscvrt.*srt2pi;
end;
if(sdepos)
%              WetDry fluxes of particles are summed
v(ityp) = dryflux + wetflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
if(sddep)
%              Dry flux of particles
v(ityp) = dryflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;
if(swdep)
%              Wet flux of particles
v(ityp) = wetflux;
if(wetscim)
vdry(ityp) = v(ityp) ./ wqcorg;
end;
end;

else;
%           Determine Deposition Correction Factors for Particles
if(~ardplete &&(ldpart || lwpart) )
[xarg, wdonly]=pdep(xarg, wdonly);
end;

%           Calculate the Vertical Term, V for particles
v(ityp) = 0.0;
if(wetscim)
vdry(ityp) = 0.0;
end;
for j = 1: npd;
%              Settling may alter SZ for the Jth particle plume
szadj = sz.*szcor(j);
a0 = -0.5./(szadj.*szadj);
%              Calculate Plume Tilt Due to Settling, HV
hv =(xarg./us) .* vgrav(j);
%              Calculate Settled Plume Height, HESETL
hesetl = he - hv;
%              Restrict settled height to be positive, so that the plume
%              does not settle below the surface -- this is the limit of
%              the tilted plume technique.
hesetl = max(0.0,hesetl);
%              Adjust Jth contribution by mass fraction and source
%              depletion
adj = phi(j) .* dqcor(j) .* wqcor(j);
if(sconc)
%                 Concentration
if(wdonly)
%                    Plume is above mixing height so set CONC = 0.0
v(ityp) = 0.0;
else;
%                    For Concentration, Complete Vertical Term is Needed for
%                    Each Particle Size Calulated at ZFLAG     ---   CALL VERT
[hesetl,szadj,a0,zflag,vj]=vert(hesetl,szadj,a0,zflag,vj);
vtmp    = adj.*pcorzr(j).*vj./szadj;
v(ityp) = v(ityp) + vtmp;
if(wetscim)
vdry(ityp) = vdry(ityp) +vtmp./wqcor(j);
end;
end;
elseif(sdepos) ;
if(wdonly)
%                    Plume is above mixing height so set DRYFLUX = 0.0
dryflux = 0.0;
else;
%                    For Dry Deposition, Complete Vertical Term is Needed for
%                    Each Particle Size Calulated at ZRDEP     ---   CALL VERT
[hesetl,szadj,a0,zrdep,vj]=vert(hesetl,szadj,a0,zrdep,vj);
%                    Calculate Dry Flux VJ/SZ
dryflux = adj.*pcorzd(j).*vdep(j).*vj./szadj;
end;
%                 Calculate Wet Flux VJ/SZ --
%                 For Wet Flux, Vertical Term is Integral of EXP terms
%                 Over All z, so VJ/SZ=SQRT(2PI)
wetflux = adj.*pscvrt(j).*srt2pi;
%                 WetDry fluxes of particles are summed
vtmp = dryflux + wetflux;
v(ityp) = v(ityp) + vtmp;
if(wetscim)
vdry(ityp) = vdry(ityp) + vtmp ./ wqcor(j);
end;
elseif(sddep) ;
if(wdonly)
%                    Plume is above mixing height so set DDEP = 0.0
v(ityp) = 0.0;
vdry(ityp) = 0.0;
else;
%                    For Dry Deposition, Complete Vertical Term is Needed for
%                    Each Particle Size Calulated at ZRDEP     ---   CALL VERT
[hesetl,szadj,a0,zrdep,vj]=vert(hesetl,szadj,a0,zrdep,vj);
%                    Calculate Dry Flux VJ/SZ
dryflux = adj.*pcorzd(j).*vdep(j).*vj./szadj;
%                    Dry flux of particles
v(ityp) = v(ityp) + dryflux;
if(wetscim)
vdry(ityp)=vdry(ityp)+dryflux./wqcor(j);
end;
end;
elseif(swdep) ;
%                 Calculate Wet Flux VJ/SZ --
%                 For Wet Flux, Vertical Term is Integral of EXP terms
%                 Over All z, so VJ/SZ=SQRT(2PI)
wetflux = adj.*pscvrt(j).*srt2pi;
%                 Wet flux of particles
v(ityp) = v(ityp) + wetflux;
if(wetscim)
vdry(ityp) = vdry(ityp)+wetflux./wqcor(j);
end;
end;
end; j = fix(npd+1);
end;

%        Calculate the Decay Term, D                        ---   CALL DECAY
[xarg]=decay(xarg);

%        Complete TERM (SZ already in denomenator of V)
rcz =(d.*v(ityp))./(srt2pi);
if(wetscim)
rczd =(d.*vdry(ityp))./(srt2pi);
end;

end;

return;
end %subroutine asimpl


function debout(varargin)
%***********************************************************************
%                 DEBOUT Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Outputs Debugging Information: Sigmas, Plume Heights,
%                 etc., for Each Calculation

%        PROGRAMMER: Roger Brode, Jeff Wang

%        DATE:    March 2, 1992

%        INPUTS:  Downwind Distance
%                 Crosswind Distance
%                 Plume Height
%                 Stack Top Wind Speed
%                 Lateral Dispersion Parameter
%                 Vertical Dispersion Parameter
%                 Stability Class
%                 Mixing Height
%                 Receptor Height Above Ground
%                 Emission Rate and Units Scaling Factor
%                 Source Parameter Arrays

%        OUTPUTS: Debug Outputs

%        CALLED FROM:   PSIMPL, AREAIN
%***********************************************************************

%     Variable Declarations
% use main1;
persistent modnam ; 

if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;

%     Variable Initializations
modnam = 'DEBOUT';

writef(fid_iounit,[ '\n ' ,repmat(' ',1,1),'JDAY= ','%3i','  IHOUR= ','%5i','  KURDAT= ','%10i','  ISRC= ','%3i','  IREC= ','%3i' ' \n'], jday, ihour, kurdat, isrc, irec);
%format[1X,'JDAY= ',i3,'  IHOUR= ',i5,'  KURDAT= ',i10,'  ISRC= ',i3,'  IREC= ',i3);

if(strcmp(deblank(srctyp(isrc)),deblank('point')))
writef(fid_iounit,[repmat(' ',1,1),'QS= ','%8.2f',' HS= ','%8.2f',' TS= ','%8.2f',' VS= ','%8.2f',' DS= ','%8.2f',' DSBH= ','%8.2f',' DSBW= ','%8.2f',' US= ','%8.5f' ' \n'], qs, hs, ts, vs, ds, dsbh, dsbw, us);
%format(1X,'QS= ',f8.2,' HS= ',f8.2,' TS= ',f8.2,' VS= ',f8.2,' DS= ',f8.2,' DSBH= ',f8.2,' DSBW= ',f8.2, ' US= ',f8.5);
writef(fid_iounit,[repmat(' ',1,1),' FB= ','%11.5f',' FM= ','%11.5f' ' \n'], fb, fm);
%format(1X,' FB= ',f11.5,' FM= ',f11.5);
elseif(strcmp(deblank(srctyp(isrc)),deblank('volume'))) ;
writef(fid_iounit,[repmat(' ',1,1),'QS= ','%8.2f','  HS= ','%8.2f','  SYINIT= ','%8.2f','  SZINIT= ','%8.2f','  US= ','%8.5f' ' \n'], qs, hs, syinit, szinit, us);
%format(1X,'QS= ',f8.2,'  HS= ',f8.2,'  SYINIT= ',f8.2,'  SZINIT= ',f8.2,'  US= ',f8.5);
elseif(strcmp(deblank(srctyp(isrc)),deblank('area'))) ;

writef(fid_iounit,[repmat(' ',1,1),'QS= ','%8.2f','  HS= ','%8.2f','  XINIT= ','%8.2f','  US= ','%8.5f','  E= ','%14.8f','  SZINIT= ','%8.2f' ' \n'], qs, hs, xinit, us, e, szinit);
%format(1X,'QS= ',f8.2,'  HS= ',f8.2,'  XINIT= ',f8.2,'  US= ',f8.5,'  E= ',g14.8, '  SZINIT= ',f8.2);

end;

writef(fid_iounit,[repmat(' ',1,1),'X= ','%12.4f','  Y= ','%12.4f','  XY= ','%11.5f','  XZ= ','%11.5f','  SY= ','%11.5f','  SZ= ','%12.5f' ' \n'], x, y, xy, xz, sy, sz);
%format(1X,'X= ',f12.4,'  Y= ',f12.4,'  XY= ',f11.5,'  XZ= ',f11.5,'  SY= ',f11.5,'  SZ= ',f12.5);
if(conc)
writef(fid_iounit,[repmat(' ',1,1),'HE= ','%11.5f','  HEMWAK= ','%11.5f','  HEFLAT= ','%11.5f','  KST= ','%2i','  TA= ','%6.1f','  ZI= ','%8.2f','  V= ','%12.6f','  D= ','%12.6f' ' \n'], he, hemwak, heflat, kst, ta, zi, v(1), d);
%format(1X,'HE= ',f11.5,'  HEMWAK= ',f11.5,'  HEFLAT= ',f11.5,'  KST= ',i2,'  TA= ',f6.1,'  ZI= ',f8.2,'  V= ',e12.6, '  D= ',e12.6);
else;
writef(fid_iounit,[repmat(' ',1,1),'HE= ','%11.5f','  HEMWAK= ','%11.5f','  HEFLAT= ','%11.5f','  KST= ','%2i','  TA= ','%6.1f','  ZI= ','%8.2f','  V= ','%12.6f','  D= ','%12.6f' ' \n'], he, hemwak, heflat, kst, ta, zi, v(1), d);
%format(1X,'HE= ',f11.5,'  HEMWAK= ',f11.5,'  HEFLAT= ',f11.5,'  KST= ',i2,'  TA= ',f6.1,'  ZI= ',f8.2,'  V= ',e12.6, '  D= ',e12.6);
end;
writef(fid_iounit,[repmat(' ',1,1),'ZLB=','%11.5f','  RINIT= ','%9.4f','  ZLY= ','%9.4f','  DA= ','%8.6f','  WAKE= ','%3f','  WAKESS=','%3f' ' \n'], zlb, rinit, zly, da, wake, wakess);
%format(1X,'ZLB=',f11.5,'  RINIT= ',f9.4,'  ZLY= ',f9.4,'  DA= ',f8.6,'  WAKE= ',l3,'  WAKESS=',l3);
writef(fid_iounit,[repmat(' ',1,1),'QTK= ','%12.5f','  XF= ','%9.2f','  XFB= ','%9.2f','  XFM= ','%9.2f','  DHF= ','%9.2f','  DHP= ','%9.2f' ' \n'], qtk, xf, xfb, xfm, dhf, dhp);
%format(1X,'QTK= ',e12.5,'  XF= ',f9.2,'  XFB= ',f9.2,'  XFM= ',f9.2,'  DHF= ',f9.2, '  DHP= ',f9.2);

writef(fid_iounit,[repmat(' ',1,1),'*** HRVAL= ','%16.8f',' ***' ' \n'], hrval(1));
%format(1X,'*** HRVAL= ',g16.8,' ***');

return;
end %subroutine debout

%----------------------------------------------------------------------

function vdp(varargin)
%----------------------------------------------------------------------

% --- ISC2ST     Version:  1.0     Level:  930215                   VDP
%                J. Scire, SRC

% --- MODIFIED   December 29, 1997
%                Removed assignment of deposition reference height, ZRDEP,
%                which is now assigned a value of 1.0m in SUB. VDP1.
%                R. W. Brode, PES, Inc.

% --- MODIFIED   May 26, 1995
%                Modified atmospheric resistance term, ra, based on
%                D. Byun and R. Dennis, Atmos. Environ., Vol. 29, No. 1
%                R. W. Brode, PES, Inc.

% --- MODIFIED   March 9, 1994
%                Changed procedure for estimating the deposition layer
%                resistance.
%                D.T. Bailey, USEPA

% --- PURPOSE:  Compute particle deposition velocities for each size
%               category of a size distribution.

% --- INPUTS:
%     Common block /METVAR/ variables:
%               Z0M - real       - Surface roughness length (m)
%             USTAR - real       - Friction velocity (m/s)
%                EL - real       - Monin-Obukhov length (m)
%     Common block /CALCS3/ variables:
%               NPD - integer    - Number of particle size categories
%             PDIAM - real array - Mean diameter (microns) of each
%                                  particle size category
%               PHI - real array - Mass fraction in each size category
%             PDENS - real       - Particle density (g/cm**3)
%                SC - real array - Schmidt number
%             VGRAV - real array - Gravitational settling velocity (m/s)
%             TSTOP - real array - Stopping time (s)
%     Common block /SOURC4/ variables:
%            VAIRMS - real       - Viscosity of air (m**2/s)
%             ZRDEP - real       - Reference height (m)
%            VDPHOR - real       - Phoretic effects term (m/s)

% --- OUTPUT:
%     Common block /CALCS3/ variables:
%              VDEP - real array - Deposition velocity (m/s) for each
%                                  particle size category

% --- VDP called by:  PCALC, VCALC, ACALC
% --- VDP calls:      none
%----------------------------------------------------------------------

% use main1;


savemlv;
integer :: :: i, io6, n;
real :: :: ustarr, ell, elabs, psih, ra, a1, b1, t1, st, xinert,schmidt, rd(npdmax), rdg, rg, b, fr, rs, rf, rc;

io6=iounit;

% ***
if(debug)
writef(fid_io6,['%s %0.15g \n'],'IHOUR  = ',ihour);
writef(fid_io6,['%s %0.15g \n'],'ISTAHR = ',istahr);
writef(fid_io6,['%s %0.15g \n'],'IENDHR = ',iendhr);
writef(fid_io6,['%s %0.15g \n'],'IEVENT = ',ievent);
writef(fid_io6,['%s %0.15g \n'],'USTAR(IHOUR) = ',austar(ihour));
writef(fid_io6,['%0.15g \n']);
writef(fid_io6,['%s \n'],'SUBR. VDP -- Inputs');
writef(fid_io6,['%s %0.15g \n'],'USTAR (m/s)     = ',ustar);
writef(fid_io6,['%s %0.15g \n'],'MONIN-EL (m)    = ',el);
writef(fid_io6,['%s %0.15g \n'],'Z0M (m)         = ',z0m);
writef(fid_io6,['%s %0.15g \n'],'VDPHOR (m/s)    = ',vdphor);
writef(fid_io6,['%s %0.15g \n'],'NPD             = ',npd);
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'PDIAM (um)      = ',pdiam(n)); end;
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'FRACT           = ',phi(n)); end;
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'PDENS (g/cm**3) = ',pdens(n)); end;
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'SC              = ',sc(n)); end;
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'VGRAV (m/s)     = ',vgrav(n)); end;
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'TSTOP (s)       = ',tstop(n)); end;
writef(fid_io6,['%s %0.15g \n'],'VAIRMS (m**2/s) = ',vairms);
writef(fid_io6,['%s %0.15g \n'],'ZRDEP (m)       = ',zrdep);
writef(fid_io6,['%s %0.15g \n'],'VDPHOR (m/s)    = ',vdphor);
end;
% ***

% --- use minimum value of USTAR to avoid numerical problems
% --- when USTAR near zero
ustarr=max(ustar,1.0e-9);

% --- Minimum absolute value of Monin-Obukhov length is 1.0 m
if(el >= 0.0)
% ---    stable
ell=max(el,1.0);
else;
% ---    unstable
ell=min(el,-1.0);
end;

% --- Calculate atmospheric resistance (s/m)
elabs=abs(ell);
if(ell > 0.0)
% ---    Stable
% ---    VK is the von Karman constant, set as parameter in MAIN1
psih = 4.7.*zrdep./ell;
ra =(1.0./(vk.*ustarr)) .*(log(zrdep./z0m) + psih);

else;
% ---    Unstable
a1 = 16..*zrdep./elabs;
b1 = 16..*z0m./elabs;
ra =(1.0./(vk.*ustarr)) .*(1.0.*log(((2.+a1)-2..*sqrt(1.+a1)) .*((2.+b1)+2..*sqrt(1.+b1)) ./(a1.*b1) ));
end;

% ***
if(debug)
writef(fid_io6,['%0.15g \n']);
writef(fid_io6,['%s %0.15g \n'],'USTARR (m/s)    = ',ustarr);
writef(fid_io6,['%s %0.15g \n'],'ELL (m)         = ',ell);
writef(fid_io6,['%s %0.15g \n'],'PSIH            = ',psih);
end;
% ***

if(npd == 0 && luservd)

% ---    GAS DEPOSITION with User-specified Deposition Velocity

vdepg = uservd;
% ***
if(debug)
writef(fid_io6,['%0.15g \n']);
writef(fid_io6,['%s \n'],'User-specified deposition velocity:');
writef(fid_io6,['%s %0.15g \n'],'VDEPG (m/s) = ',vdepg);
end;
% ***

elseif(npd == 0) ;

% ---    GAS DEPOSITION

% ---    Compute the deposition layer resistance for gases, RDG
% ---    (RD1 is d1*SC**d2/vk -- computed in setup routine)
rdg=rd1(isrc)./ustarr;

% ---    Compute resistance directly to ground or water, RG
% ---    Water is assumed if LAI = 0.0
if(xlai == 0.0)

% ---       Water cell (RGW1 computed in setup routine as
% ---                   RGW1 = HENRY/(ALPHAS * D3)
rg=rgw1(isrc)./ustarr;
else;

% ---       Land cell (RG computed in setup routine)
rg=rgg(isrc);
end;

% ---    Stomatal pore resistance (RS)
if(unstressed)

% ---       Vegetation is activeunstressed (IVEG=1 in CALPUFF)
%           (B = stomatal pore opening (m), BMIN = minimum stomatal
%           opening, BMAX = maximum stomatal opening, fr is the approx.
%           fraction of peak short-wave solar radiation available for a
%           particular hour)

% ---       Temperature effects -- If T > 35 deg. C, stomata fully
% ---       open to allow evaporative cooling -- but only if unstressed --
% ---       If T < 10 deg. C, stomata closed due to decreased metabolic
% ---       activity)
if(ta > 308.)
% ---          T > 35 deg. C
b=bmax;
elseif(ta < 283.);
% ---          T < 10 deg. C
b=bmin;
else;
fr=qsw./qswmax;
fr=max(0.0,fr);
fr=min(1.0,fr);
b=bmax.*fr+bmin.*(1.-fr);
end;

rs=pconst./(b.*pdiff(isrc));
elseif(stressed) ;

% ---       Vegetation is active and stressed (IVEG=2 in CALPUFF)
%           (Stomatal opening is at its minimum size)
rs=pconst./(bmin.*pdiff(isrc));
elseif(inactive) ;

% ---       Vegetation is inactive (IVEG=3)
rs=9.9e9;
end;

% ---    Internal foliage resistance (RF)
%        (RM is the mesophyll resistance)
rf=rs+rm(isrc);

% ---    Compute canopy resistance
rc=1.0./(xlai./rf+xlai./rcut(isrc)+1.0./rg);

% ---    Deposition velocity is the inverse of the sum of the
% ---    atmospheric, deposition layer, and canopy resistances
vdepg = 1.0./(ra+rdg+rc);

% ***
if(debug)
writef(fid_io6,['%0.15g \n']);
writef(fid_io6,['%s %0.15g \n'],'RA (s/m)    = ',ra);
writef(fid_io6,['%s %0.15g \n'],'RDG (s/m)   = ',rdg);
writef(fid_io6,['%s %0.15g \n'],'RC (s/m)    = ',rc);
writef(fid_io6,['%s %0.15g \n'],'VDEPG (m/s) = ',vdepg);
end;
% ***

else;


% ---    PARTICLE DEPOSITION

t1=ustarr.*ustarr./vairms;

% ---    LOOP OVER SIZE INTERVALS
for i=1:npd;

st=tstop(i).*t1;

% ---       Compute inertial impaction term
xinert=10.^(-3../st);

% ---       Adjust (raise) the Schmidt Number to the 2/3rd's power.
schmidt = sc(i) .^(-.667)                                       dtb94068;

% ---       Compute the deposition layer resistance (s/m)
rd(i)=1.0 ./(ustarr .*(schmidt + xinert))                        dtb94068;

% ---       Deposition velocity for this current interval
vdep(i)=1.0./(ra+rd(i)+ra.*rd(i).*vgrav(i))+vgrav(i)+vdphor;

end; i=npd+1;
% ***
if(debug)
writef(fid_io6,['%0.15g \n']);
writef(fid_io6,['%s %0.15g \n'],'RA (s/m)    = ',ra);
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'RD (s/m)    = ',rd(n)); end;
for n=(1):(npd), writef(fid_io6,['%s %0.15g \n'],'VDEP (m/s)  = ',vdep(n)); end;
end;
% ***

end;

return;
end %subroutine vdp

%-----------------------------------------------------------------------

function setszmn(varargin)
%-----------------------------------------------------------------------

% --- ISCST2    Version: 1.0            Level: 931215           SETSZMN
%               D. Strimaitis, SRC

% PURPOSE:     SETSZMN determines the value of sigma-z at which the rate
%              of growth in sigma-z equals the rate at which the settling
%              velocity acts to reduce the height of the center-of-mass.
%              A default minimum of 2*zd, where zd is the near-surface
%              height at which the deposition flux is evaluated, is
%              returned if there is no balance-point.

% ARGUMENTS:  (MAIN1)
%    PASSED:  kst       stability class (A=1, F=6)                   [i]
%             zrdep     reference height for deposition flux  (m)    [r]
%             vs        settling velocity  (m/s)                     [r]
%             us        plume advection wind speed (m/s)             [r]
%             urban     logical for URBAN/RURAL dispersion params    [l]
%             npd       number of particle size categories           [i]

%  RETURNED:  szmin     Minimum value of sigma-z (m)                 [r]

% CALLING ROUTINES:   PCALC, VCALC, ACALC

% EXTERNAL ROUTINES:  GCUBIC
%-----------------------------------------------------------------------
% use main1;

persistent car cau cbr cbu firstCall root ; if isempty(firstCall),firstCall=1;end; 
a1=[];a2=[];a3=[];
savemlv;
integer :: :: i, j;
real :: :: xmin, a, b, c, a1, a2, a3, aby2csq;

if isempty(root), root=zeros(1,3); end;
if isempty(car), car=zeros(1,6); end;
if isempty(cau), cau=zeros(1,6); end;
if isempty(cbr), cbr=zeros(1,6); end;
if isempty(cbu), cbu=zeros(1,6); end;

if firstCall,   car=[.2,.12,.08,.06,.03,.016];  end;
if firstCall,   cbr=[0.,0.,.0002,.0015,.0003,.0003];  end;
if firstCall,   cau=[.24,.24,.2,.14,.08,.08];  end;
if firstCall,   cbu=[.001,.001,0.,.0003,.0015,.0015];  end;
firstCall=0;


% --- Loop over particle sizes
for i=1:npd;
xmin=0.0;
szmin(i)=2..*zrdep;
c=rtpiby2.*vgrav(i)./us;

% ---    Urban section
if(urban)
a=cau(kst);
b=cbu(kst);
if(kst >= 4)
if(a > 20..*c)
szmin(i)=a.*a./(2..*b.*c);
elseif(a > c) ;
% ---             Solve cubic for y=bx, then report x       ---  call GCUBIC
aby2csq=(a./(2..*c)).^2;
a1=(3.-aby2csq);
a2=(3.-4..*aby2csq);
a3=(1.-4..*aby2csq);
[a1,a2,a3,root]=gcubic(a1,a2,a3,root);
% ---             There should be ONE real root
if(root(2) ~= 0. || root(3) ~= 0.)
writef(1,['%s \n'], 'SETSZMN: Potential error!!! ');
writef(1,['%s \n'], 'More than 1 root ----');
for j=(1):(3), writef(1,['%s %0.15g \n'], 'xb= ',root(j)); end;
end;
xmin=root(1)./b;
szmin(i)=a.*xmin./sqrt(1.+b.*xmin);
end;
end;

% ---    Rural section
else;
a=car(kst);
b=cbr(kst);
if(kst == 3 || kst == 4)
if(a > 20..*c)
szmin(i)=a.*a./(2..*b.*c);
elseif(a > c) ;
% ---             Solve cubic for y=bx, then report x       ---  call GCUBIC
aby2csq=(a./(2..*c)).^2;
a1=(3.-aby2csq);
a2=(3.-4..*aby2csq);
a3=(1.-4..*aby2csq);
[a1,a2,a3,root]=gcubic(a1,a2,a3,root);
% ---             There should be ONE real root
if(root(2) ~= 0. || root(3) ~= 0.)
writef(1,['%s \n'], 'Potential error!!! More than 1 root');
for j=(1):(3), writef(1,['%s %0.15g \n'], 'xb= ',root(j)); end;
end;
xmin=root(1)./b;
szmin(i)=a.*xmin./sqrt(1.+b.*xmin);
end;
elseif(kst > 4) ;
if(a > c)
xmin=(sqrt(a./c)-1.)./b;
szmin(i)=a.*xmin./(1+b.*xmin);
end;
end;
end;

end; i=npd+1;

return;
end %subroutine setszmn

%-----------------------------------------------------------------------

function [a1,a2,a3,root]=gcubic(a1,a2,a3,root);
%-----------------------------------------------------------------------

% --- ISCST2    Version: 1.0            Level: 931215           GCUBIC
%               D. Strimaitis, SRC

% PURPOSE:     Program solves the general cubic equation of the form:
%                  0 = x**3 + (a1)x**2 + (a2)x + (a3)
%              for the real roots
%              (Numerical Recipes, Press et al., 1986)

% ARGUMENTS:
%    PASSED:  a1,a2,a3  constants for terms as described above       [r]

%  RETURNED:  root      root(s) of equation                          [r]

% CALLING ROUTINES:   (utility routine)

% EXTERNAL ROUTINES:  none
%-----------------------------------------------------------------------



no type, intent(in)                      :: a1;

real :: :: a1, third, a1sq, a1cube, a1by3, q, r, qcube,rsq, sqrtq2, theta, arg;

real,:: parameter :: twopi=6.2831853, fourpi=12.566371;

third=1../3.;
a1sq=a1.*a1;
a1cube=a1.*a1sq;
a1by3=a1.*third;

q=(a1sq-3..*a2)./9.;
r=(2..*a1cube-9..*a1.*a2+27..*a3)./54.;

qcube=q.*q.*q;
rsq=r.*r;

if(qcube >= rsq)
% ---    THREE real roots
sqrtq2=sqrt(q).*2.;
theta=acos(r./sqrt(qcube));
root(1)=-sqrtq2.*cos(theta./3.)-a1by3;
root(2)=-sqrtq2.*cos((theta+twopi)./3.)-a1by3;
root(3)=-sqrtq2.*cos((theta+fourpi)./3.)-a1by3;
else;
% ---    ONE real root
arg=(sqrt(rsq-qcube)+abs(r)).^third;
root(1)=-(abs(1.0).*sign(r)).*(arg+q./arg)-a1by3;
root(2)=0.;
root(3)=0.;
end;


return;
end %subroutine gcubic

%----------------------------------------------------------------------

function scavrat(varargin)
%----------------------------------------------------------------------

% --- ISCST2     Version: 1.0       Level: 931108               SCAVRAT
%                D. Strimaitis, SRC

% --- PURPOSE:  Compute the wet SCAVenging RATio for particles, as a
%               function of particle size, and for gases

% --- INPUTS:
%     Common block /METVAR/ variables:
%            IPCODE - integer    - Precip. code (00-45)
%             PRATE - real       - Precip. rate (mm/hr)
%                TA - real       - Ambient Temperature (deg K)
%     Common block /CALCS3/ variables:
%               NPD - integer    - Number of particle size categories
%             PSCAV - real array - Particle scavenging coefs. for liquid
%                                  (1) and frozen (2) precip. for each
%                                  size category (1/[s-mm/hr])
%             GSCAV - real array - Gas scavenging coefs. for liquid (1)
%                                  and frozen (2) precip. (1/[s-mm/hr])

% --- OUTPUT:
%     Common block /CALCS3/ variables:
%            PSCVRT - real array - Scavenging ratio for particles (1/s)
%            GSCVRT - real       - Scavenging ratio for gases (1/s)

% --- SCAVRAT called by:  PCALC, VCALC, ACALC
% --- SCAVRAT calls:      none
%----------------------------------------------------------------------

% --- Include common blocks
% use main1;

persistent firstCall imiss ; if isempty(firstCall),firstCall=1;end; 

savemlv;
integer :: :: i, n, ilq, imiss;

if firstCall,   imiss=[9999];  end;
firstCall=0;

if(debug)
writef(fid_iounit,['%0.15g \n']);
writef(fid_iounit,['%s \n'],'SUBR. SCAVRAT -- Inputs');
writef(fid_iounit,['%s %0.15g \n'],'IPCODE               = ',ipcode);
writef(fid_iounit,['%s %0.15g \n'],'PRATE (mm/hr)        = ',prate);
writef(fid_iounit,['%s %0.15g \n'],'TA (deg K)           = ',ta);
writef(fid_iounit,['%s %0.15g \n'],'NPD                  = ',npd);
for n=(1):(npd), writef(fid_iounit,['%s %0.15g \n'],'PSCAV(1) 1/(s-mm/hr) = ',pscav(n,1)); end;
for n=(1):(npd), writef(fid_iounit,['%s %0.15g \n'],'PSCAV(2) 1/(s-mm/hr) = ',pscav(n,2)); end;
writef(fid_iounit,['%s %0.15g \n'],'GSCAV(1) 1/(s-mm/hr) = ',gscav(1));
writef(fid_iounit,['%s %0.15g \n'],'GSCAV(2) 1/(s-mm/hr) = ',gscav(2));
writef(fid_iounit,['%s \n'],' (1 = Liquid ; 2 = Frozen )');
writef(fid_iounit,['%0.15g \n']);
end;

% --- If no precipitation, no wet removal
if(prate == 0.)
for i=1:npd;
pscvrt(i)=0.0;
end; i=npd+1;
gscvrt=0.0;
else;
% ---    Determine if precip. is liquid (ILQ=1) or frozen (ILQ=2)
if(ipcode == imiss || ipcode == 0)
% ---       Precip. code is unavailable due to missing data or no
%           precip. at time of obs. at surface station, therefore,
%           determine precip. type based on the air temperature
% ---       Assume liquid precip. if temp. > freezing, otherwise,
%           assume frozen precip.
if(ta > 273.15)
ilq=1;
else;
ilq=2;
end;
elseif(ipcode <= 18) ;
% ---       Liquid precipitation type
ilq=1;
else;
% ---       Frozen precipitation type
ilq=2;
end;
% ---    Determine the scavenging ratios
for i=1:npd;
pscvrt(i)=pscav(i,ilq).*prate;
end; i=npd+1;
gscvrt=gscav(ilq).*prate;
end;

if(debug)
writef(fid_iounit,['%s \n'],'SUBR. SCAVRAT -- Results');
writef(fid_iounit,['%s %0.15g \n'],'GSCVRT (1/s)= ',gscvrt);
for n=(1):(npd), writef(fid_iounit,['%s %0.15g \n'],'PSCVRT (1/s)= ',pscvrt(n)); end;
writef(fid_iounit,['%0.15g \n']);
end;

return;
end %subroutine scavrat


function [xarg, lwdonly]=pdep(xarg, lwdonly);
%***********************************************************************
%               PDEP Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Simple Terrain Deposition Adjustment
%                 Factors from DEPCOR

%        PROGRAMMER: R. W. Brode, PES, Inc.

%        DATE:    September 30, 1994

%        MODIFIED:   To add logical argument for Wet Deposition Only,
%                    to skip call to DEPCOR when plume is above ZI.
%                    R.W. Brode, PES, 7/17/95

%        INPUTS:     LWDONLY, logical specifying whether Wet Deposition
%                    Only is to be calculated for plume above ZI

%        OUTPUTS:


%        CALLED FROM:   PSIMPL
%***********************************************************************

%     Variable Declarations
% use main1;

persistent i lterr modnam ; 
vdep=[];vgrav=[];zrdep=[];zflag=[];xz=[];heflat=[];zi=[];us=[];xs=[];ys=[];xr=[];yr=[];rural=[];urban=[];kst=[];sz=[];sbid=[];szmin=[];zelev=[];zs=[];debug=[];iounit=[];srctyp=[];isrc=[];ltgrid=[];kurdat=[];dqcor=[];pcorzr=[];pcorzd=[];szcor=[];toxics=[];
no type, intent(in)                      :: xarg;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(i), i=0; end;
if isempty(lterr), lterr=false; end;

%     Variable Initializations
modnam = 'PDEP';

%     Set LTERR to falsemlv to signal simple terrain call to DEPCOR.
lterr = false;

%     Loop over particle sizes
for i=1:npd;
dqcor(i)  = 1.0;
pcorzr(i) = 1.0;
pcorzd(i) = 1.0;
szcor(i)  = 1.0;
%        Initialize wetdry source depletion factors,
%        profile correction factors, and settles sigma-z
%        factors to unity. - Done in DEPCOR
if(ddplete && ~lwdonly)
%           Determine factors for depletion - note that
%           plume ht adjustment for terrain is signalled
%           by a local logical - LTERR
%           Simple Terrain Model          ---   CALL DEPCOR
[ vdep(i),vgrav(i),zrdep,zflag,xarg,xz,heflat,zi,us,xs,ys,xr,yr, rural,urban,kst,sz,sbid,szmin(i),zelev,zs,lterr,debug,iounit, srctyp(isrc),ltgrid,kurdat,dqcor(i),pcorzr(i),pcorzd(i),szcor(i),toxics]=depcor( vdep(i),vgrav(i),zrdep,zflag,xarg,xz,heflat,zi,us,xs,ys,xr,yr, rural,urban,kst,sz,sbid,szmin(i),zelev,zs,lterr,debug,iounit, srctyp(isrc),ltgrid,kurdat,dqcor(i),pcorzr(i),pcorzd(i),szcor(i),toxics);
end;
if(wdplete)
%           Determine source depletion factor
%           from wet removal
%           Simple Terrain Model
wqcor(i) = exp(-pscvrt(i).*xarg./us);
else;
wqcor(i) = 1.;
end;
end; i=fix(npd+1);

return;
end %subroutine pdep


function [lwdonly]=pdepc(lwdonly);
%***********************************************************************
%               PDEPC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Complex Terrain Deposition Adjustment
%                 Factors from DEPCOR

%        PROGRAMMER: R. W. Brode, PES, Inc.

%        DATE:    September 30, 1994

%        MODIFIED:   To add logical argument for Wet Deposition Only,
%                    to skip call to DEPCOR when plume is above ZI.
%                    R.W. Brode, PES, 7/17/95

%        INPUTS:     LWDONLY, logical specifying whether Wet Deposition
%                    Only is to be calculated for plume above ZI

%        OUTPUTS:


%        CALLED FROM:   PCOMPL
%***********************************************************************

%     Variable Declarations
% use main1;

persistent i lterr modnam ; 
vdep=[];vgrav=[];zrdep=[];zflag=[];distr=[];xzcmp1=[];hecomp=[];zi=[];us=[];xs=[];ys=[];xr=[];yr=[];rural=[];urban=[];kst=[];szcmp1=[];sbcmp1=[];szmin=[];zelev=[];zs=[];debug=[];iounit=[];srctyp=[];isrc=[];ltgrid=[];kurdat=[];dqcorc=[];pcorzrc=[];pcorzdc=[];szcorc=[];toxics=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(i), i=0; end;
if isempty(lterr), lterr=false; end;

%     Variable Initializations
modnam = 'PDEPC';

%     Set LTERR to truemlv to signal complex terrain call to DEPCOR.
lterr = true;

%     Loop over particle sizes
for i=1:npd;
dqcorc(i)  = 1.0;
pcorzrc(i) = 1.0;
pcorzdc(i) = 1.0;
szcorc(i)  = 1.0;
%        Initialize wetdry source depletion factors,
%        profile correction factors, and settles sigma-z
%        factors to unity. - Done in DEPCOR
if(ddplete && ~lwdonly)
%           Determine factors for depletion - note that
%           plume ht adjustment for terrain is signalled
%           by a local logical - LTERR
%           Complex Terrain Model         ---   CALL DEPCOR
[ vdep(i),vgrav(i),zrdep,zflag,distr,xzcmp1,hecomp,zi,us,xs,ys,xr,yr, rural,urban,kst,szcmp1,sbcmp1,szmin(i),zelev,zs,lterr,debug,iounit, srctyp(isrc),ltgrid,kurdat,dqcorc(i),pcorzrc(i),pcorzdc(i), szcorc(i),toxics]=depcor( vdep(i),vgrav(i),zrdep,zflag,distr,xzcmp1,hecomp,zi,us,xs,ys,xr,yr, rural,urban,kst,szcmp1,sbcmp1,szmin(i),zelev,zs,lterr,debug,iounit, srctyp(isrc),ltgrid,kurdat,dqcorc(i),pcorzrc(i),pcorzdc(i), szcorc(i),toxics);
end;
if(wdplete)
%           Determine source depletion factor
%           from wet removal
%           Complex Terrain Model - use radial distance
wqcorc(i) = exp(-pscvrt(i).*distr./us);
else;
wqcorc(i) = 1.;
end;
end; i=fix(npd+1);

return;
end %subroutine pdepc

function [xarg, lwdonly]=pdepg(xarg, lwdonly);
%***********************************************************************
%               PDEPG Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Simple Terrain Deposition Adjustment
%                 Factors from DEPCOR for Gases

%        PROGRAMMER: R. W. Brode, PES, Inc.

%        DATE:       May 6, 1996

%        INPUTS:     LWDONLY, logical specifying whether Wet Deposition
%                    Only is to be calculated for plume above ZI

%        OUTPUTS:


%        CALLED FROM:   PSIMPL
%***********************************************************************

%     Variable Declarations
% use main1;

persistent lterr modnam ; 
vdepg=[];zrdep=[];zflag=[];xz=[];heflat=[];zi=[];us=[];xs=[];ys=[];xr=[];yr=[];rural=[];urban=[];kst=[];sz=[];sbid=[];zelev=[];zs=[];debug=[];iounit=[];srctyp=[];isrc=[];ltgrid=[];kurdat=[];dqcorg=[];pcorzrg=[];pcorzdg=[];szcorg=[];toxics=[];
no type, intent(in)                      :: xarg;
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(lterr), lterr=false; end;

%     Variable Initializations
modnam = 'PDEPG';

%     Set LTERR to falsemlv to signal simple terrain call to DEPCOR.
lterr = false;

%     Initialize source depletion factors to unity.
dqcorg  = 1.0;
pcorzrg = 1.0;
pcorzdg = 1.0;
szcorg  = 1.0;
wqcorg  = 1.0;
%     Initialize wetdry source depletion factors,
%     profile correction factors, and settles sigma-z
%     factors to unity. - Done in DEPCOR
if(ddplete && ~lwdonly)
%        Determine factors for depletion - note that
%        plume ht adjustment for terrain is signalled
%        by a local logical - LTERR
%        Simple Terrain Model          ---   CALL DEPCOR
[ vdepg,dumvar2,zrdep,zflag, xarg,xz,heflat,zi,us,xs,ys,xr,yr,rural,urban,kst,sz,sbid,dumvar19,zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat, dqcorg,pcorzrg,pcorzdg,szcorg,toxics]=depcor( vdepg,0.0,zrdep,zflag, xarg,xz,heflat,zi,us,xs,ys,xr,yr,rural,urban,kst,sz,sbid, 2..*zrdep,zelev,zs,lterr,debug,iounit,srctyp(isrc),ltgrid,kurdat, dqcorg,pcorzrg,pcorzdg,szcorg,toxics);
end;
if(wdplete)
%        Determine source depletion factor
%        from wet removal (GASES)
%        Simple Terrain Model
wqcorg = exp(-gscvrt.*xarg./us);
end;

return;
end %subroutine pdepg

function [lwdonly]=pdepgc(lwdonly);
%***********************************************************************
%               PDEPGC Module of ISC2 Short Term Model - ISCST2

%        PURPOSE: Calculates Complex Terrain Deposition Adjustment
%                 Factors from DEPCOR for Gases

%        PROGRAMMER: R. W. Brode, PES, Inc.

%        DATE:       May 6, 1996

%        INPUTS:     LWDONLY, logical specifying whether Wet Deposition
%                    Only is to be calculated for plume above ZI

%        OUTPUTS:


%        CALLED FROM:   PCOMPL
%***********************************************************************

%     Variable Declarations
% use main1;

persistent lterr modnam ; 
vdepg=[];zrdep=[];zflag=[];distr=[];xzcmp1=[];hecomp=[];zi=[];us=[];xs=[];ys=[];xr=[];yr=[];rural=[];urban=[];kst=[];szcmp1=[];sbcmp1=[];zelev=[];zs=[];debug=[];iounit=[];srctyp=[];isrc=[];ltgrid=[];kurdat=[];dqcorgc=[];pcorzrgc=[];pcorzdgc=[];szcorgc=[];toxics=[];
if isempty(modnam), modnam=repmat(' ',1,12); end;

savemlv;
if isempty(lterr), lterr=false; end;

%     Variable Initializations
modnam = 'PDEPGC';

%     Set LTERR to truemlv to signal complex terrain call to DEPCOR.
lterr = true;

%     Initialize source depletion factors to unity.
dqcorgc  = 1.0;
pcorzrgc = 1.0;
pcorzdgc = 1.0;
szcorgc  = 1.0;
wqcorgc  = 1.0;
%     Initialize wetdry source depletion factors,
%     profile correction factors, and settles sigma-z
%     factors to unity. - Done in DEPCOR
if(ddplete && ~lwdonly)
%        Determine factors for depletion - note that
%        plume ht adjustment for terrain is signalled
%        by a local logical - LTERR
%        Simple Terrain Model          ---   CALL DEPCOR
[ vdepg,dumvar2,zrdep,zflag,distr,xzcmp1,hecomp,zi,us,xs,ys,xr,yr, rural,urban,kst,szcmp1,sbcmp1,dumvar19,zelev,zs,lterr,debug,iounit, srctyp(isrc),ltgrid,kurdat,dqcorgc,pcorzrgc,pcorzdgc,szcorgc,toxics]=depcor( vdepg,0.0,zrdep,zflag,distr,xzcmp1,hecomp,zi,us,xs,ys,xr,yr, rural,urban,kst,szcmp1,sbcmp1,2..*zrdep,zelev,zs,lterr,debug,iounit, srctyp(isrc),ltgrid,kurdat,dqcorgc,pcorzrgc,pcorzdgc,szcorgc,toxics);
end;
if(wdplete)
%        Determine source depletion factor
%        from wet removal (GASES)
%        Simple Terrain Model
wqcorgc = exp(-gscvrt.*distr./us);
end;

return;
end %subroutine pdepgc




function out=writef(fid,varargin)
% function out=writef(fid,varargin)
%  Catches fortran stdout (6) and reroutes in to Matlab's stdout (1)
%  Catches fortran stderr (0) and reroutes in to Matlab's stderr (2)
if isnumeric(fid)
 if fid==6,      out=fprintf(1,varargin{:});
 elseif fid==0,  out=fprintf(2,varargin{:});
 elseif isempty(fid) %% treat empty array like a string array [sethg 2008-03-03]
  out=sprintf(varargin{:});
  if nargin>2 %set the calling var to out
   if ~isempty(inputname(1)), assignin('caller',inputname(1),out); end
  end
 else,           out=fprintf(fid,varargin{:});
 end
elseif ischar(fid)
 out=sprintf(varargin{:});
 if nargin>2 %set the calling var to out
  if ~isempty(inputname(1)), assignin('caller',inputname(1),out); end
 end
else,            out=fprintf(fid,varargin{:});
end
end