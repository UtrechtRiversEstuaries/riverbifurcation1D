%1D morphodynamical model with gradually varied flow
%for bifurcated and confluencing river networks
%by Maarten Kleinhans, April 2009
%provides input for takke43.m model


%%
%switches in model
switchC      = 2;	%if 1 then constant C, if 2 then White-Colebrook_rough constant kc
switchtransp = 1;   %if 1 then Engelund-Hansen, if 2 then MPM_adapted
switchnodal  = 3;   %if 1 then Wang, if 2 then Bolla, if 3 then MEANDERING
switchBvar   = 1;   %if 1 then constant, if 2 then variable downstream widths hydraulic geometry
switchRvar   = 3;   %if 1 then constant bend radius, 2 switching meander migration, 3 sinusoidal migration
switchCour   = 2;   %if 1 then constant time step, if 2 then critical courant adaptation

%switches for figures (1=make, 0=do not make)
switchfig1 = 1; 	%time series of water and bed levels
switchfig2 = 1; 	%time series of relative discharge in bifurcated branches
switchfig3 = 0; 	%comparison of water depths
switchfig4 = 1; 	%-year evolution of long profile
switchfig5 = 0; 	%nodal point representation Qs2/Qs3 vs Q2/Q3
switchfig6 = 0; 	%time-series downstream widths


%%
%INPUT sizes and topology

%upstream connections: (none if upstream boundary, one if bifur, two if confluence)
%downstream connections: (none if sea, one if confluence, two if bifurcation)
%confluences: upstream branches MUST have same order, if necessary add branch

%Orde = from.. to.. branches in same bifurcation 'order'
Orde = [...
    1 1;...1
    2 3;...2
    ];

%Topo = (upstr1,upstr2,downstr1,downstr2) for each branch
Topo = [...
    NaN NaN 2   3   ; ...1
    1   NaN NaN NaN ; ...2
    1   NaN NaN NaN ; ...3
    ];

%Sizes = (Q,L,D50 (mm),Mudpercentage) for each branch
Sizes = [...
    2500   6000  2 0; ...1
    1250   6000  2 0; ...2
    1250   6000  2 0; ...3
    ];

%Activity = (begin of activity (cal yr BP),end of activity) for each branch 
Activity = [...
    25   0; ...1
    25   0; ...2
    25   0; ...3
    ];

%Bifurcations = (heightaboveNAP,Rfac,Tbend,Lbend,Abend)
%for the upstream boundary, each bifurcation
%Upstream bend: relative bend radius R/W at each bifurcation
%bifurcate 2 is outer bend for positive (for variable R switchRvar=2)
%bifurcate 3 is inner bend for positive (for variable R switchRvar=2)
%period (yr) at which bend radius changes sign, REAL period is 2Tbend!
							%(for variable R switchRvar=2 and 3 resp)
Bifurcations = [...
      1.2  NaN  NaN  NaN   NaN  ; ...upstream boundary
      0.6  1  20  8000  2000 ; ...bifurcation 1 (most upstream)
    ];

%Confluences = heightaboveNAP for confluences
Confluences = [...
      ...confluence 1 (most upstream)
    ];

%Connections = heightaboveNAP for throughflow connection nodes
Connections = [...
      ...node 1 (most upstream)
    ];

Bup = 504;        %upstream width
%Basin (local datum)
xi0 = 0; %initial downstream water level
baselevelrise = 0.000; %rate in m per year


%%
%schematisation
%space
dx = 500;			%m spatial step
au = 1.0;			%1 = full upwind, 0.5 = central difference
waterlevelprecision = 1e-3; %1e-4
dischargeprecision = 1e-1; %1e-2
maxiter = 10;       %20 maximum number of backwater iterations at a bifurcation
%Hmin = 0.8;         %minimum water depth

%time
Durat = 25;          %50 years duration
dt = 0.01;       	%initial 0.01 year time step
CourantCrit = 0.8;	%au-0.1;% critical courant number for decreasing the time step
dtmin = dt/3;		%minimum time step regardless of Courant criterion
maxNt = 3; 			%maximum exceedance factor of expected time steps

%reports
Nreport = 10;			%nr of reports for bed and water level curves (>6 cyclic colors)

%Tectonics in both bifurcates gradual uplift (+) or downwarp (-) of bed
%etatect = [0 0 0 0 0 0 0].*1e-3; %tectonics rate in m per year
%ftect = [0 0 0 0 0 0 0 ; 1 1 1 1 1 1 1]; %fraction of L for up/downstream end of 'tectonic block'
%location in xcoord for report full time series of base, water, bed levels and Cr
%frep = [0.95 0.05 0.05 0.05 0.05 0.05 0.05]; %fractions of Li for time report
%nu even automatisch goed:
etatect = [ zeros(1,length(Topo(:,1))) ];
ftect   = [ zeros(1,length(Topo(:,1))) ; ones(1,length(Topo(:,1))) ];
frep    = [ ones(1,length(Topo(:,1))) ] .*0.05;

%%
%constants
Qbreak = 1000;	%discharge at which one branch is considered closed
If = 1;				%Intermittency
lamp = 0.3;			%Bed Porosity
kc = 0.15; 		    %m Roughness Height, not necessary to calibrate for WhiteColebrook 
bB = 1;             %0.6 power on discharge Q for hydraulic geometry
bD = -0.;           %-0.25 power on diameter D for hydraulic geometry
bM = -0.;           %-0.37 power on mud percentage M for hydraulic geometry

%tons/a Imposed annual sediment transport rate fed in from upstream (which must all be carried during floods)
%Gtf = 1.625E+06;	%1e6 NOT USED

%nodal point constants
k = 1;				%2 POWER OF Wang NODAL POINT RELATION
alb = 1;			%1, r parameter in Ikeda/Bolla NODAL POINT RELATION
%calibrated on Delft3D:
alw = 2;			%2, alpha parameter in Bolla NODAL POINT RELATION
epsilon = 2;		%2, spiral flow calibration parameter
