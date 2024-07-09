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

