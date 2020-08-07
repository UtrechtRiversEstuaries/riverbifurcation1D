%1D morphodynamical model with gradually varied flow
%for bifurcated and confluencing river networks
%by Maarten Kleinhans, November 2013
%provides input for takke44.m model


%% switches in model
switchC = 1;		%1 if 1 then constant C, if 2 then White-Colebrook_rough constant kc
switchCp = 2;       %2 if 1 then total shear stress, if 2 grain-related
switchBvar = 2;     %1 if 1 then constant, if 2 then variable downstream widths hydraulic geometry
switchRvar = 3;     %1 if 1 then constant bend radius, 2 switching meander migration, 3 sinusoidal migration
switchCour = 2;     %2 if 1 then constant time step, if 2 then critical courant adaptation
switchtransp = 1;   %3 if 1 then Engelund-Hansen, if 2 then MPM_adapted, if 3 then van Rijn 2007
switchvRsusp = 0;   %1 if van Rijn this switches suspended transport on
switchnodal = 3;    %2 if 1 then Wang, if 2 then Bolla, if 3 then MEANDERING
switchtransv = 0;   %1 if 1 then use transverse bed slope r of Bolla for meandering
switchperturb = 0;  %1 if 1 then perturbation put in given location after given time

%switches for figures (1=make, 0=do not make)
switchfig1 = 0; 	%time series of water and bed levels
switchfig2 = 1; 	%time series of relative discharge in bifurcated branches
switchfig3 = 0; 	%comparison of water depths
switchfig4 = 0; 	%-year evolution of long profile
switchfig5 = 0; 	%nodal point representation Qs2/Qs3 vs Q2/Q3
switchfig6 = 1; 	%time-series downstream widths


%% INPUT sizes and topology

%upstream connections: (none if upstream boundary, one if bifur, two if confluence)
%downstream connections: (none if sea, one if confluence, two if bifurcation)
%confluences: upstream branches MUST have same order, if necessary add branch

%Orde = from.. to.. branches in same bifurcation 'order'
Orde = [...
    1 1;...1
    2 3;...2
    ];

Names = {'Ganga','Hooghly','Ganga'};

%Activity = (begin of activity (cal yr BP),end of activity) for each branch 
Activity = [...
    900  0 ; ...1 eerste
    900  0 ; ...2 tweede
    900  0 ; ...3 derde
    ];

%Topo = (upstr1,upstr2,downstr1,downstr2) for each branch
Topo = [...
    NaN NaN 2   3   ; ...1
    1   NaN NaN NaN ; ...2
    1   NaN NaN NaN ; ...3
    ];

%Sizes = (Q,L,D50 (mm),Mudpercentage) for each branch
Sizes = [...
    45000   100000  0.14  0; ...1
    40000   300000  0.14  0; ...2
     5000   300000  0.14  0; ...3
    ];
%    25000   330000  0.2 0; ...2
%     5000   390000  0.2 0; ...3

%% Elevation data
%when chckelev43 is used for fitting to power profile, 
%FOR POWER PROFILE power nS=1 (linear) or higher (adjust in chckelev43)
%take elevations of LOWEST branch and FILL IN for missing paths

%Bifurcations = (heightaboveNAP,Rfac,Tbend,Lbend,Abend)
%for the upstream boundary, each bifurcation
%Upstream bend: relative bend radius R/W at each bifurcation
%bifurcate 2 is outer bend for positive (for variable R switchRvar=2)
%bifurcate 3 is inner bend for positive (for variable R switchRvar=2)
%period (yr) at which bend radius changes sign, REAL period is 2Tbend!
							%(for variable R switchRvar=2 and 3 resp)

%migration rate 0.483 km/yr (1990-2002), meander bend length 29 km and amplitude 10 km
Bifurcations = [...
     20  NaN  NaN        NaN    NaN  ; ...upstream boundary S=.00005=5x10^-5
     15  -3   29000/483  29000  -10000 ; ...bifurcation 1 (most upstream)
    ];

%Confluences = heightaboveNAP for confluences
Confluences = [...
      ...confluence 1 (most upstream)
    ];

%Connections = heightaboveNAP for throughflow connection nodes
Connections = [...
      ...node 1 (most upstream)
    ];

Bup = 1800;        %upstream width
%Basin (local datum)
xi0 = 0; %initial downstream water level
baselevelrise = 0.000; %rate in m per year


%% schematisation
%space
dx = 2500;			%m spatial step
au = 1.0;			%1 = full upwind, 0.5 = central difference
waterlevelprecision = 1e-4; %1e-4
dischargeprecision = 1e-2; %1e-2
maxiter = 10;       %20 maximum number of backwater iterations at a bifurcation
%netiter = 10;        %1 default number of network iterations netiter*Norder

%time
Durat = 300;		%900 max years duration
dt = 0.01;			%initial 0.01 year time step
CourantCrit = 0.7;	%au-0.1;% critical courant number for decreasing the time step
dtmin = dt/3;		%minimum time step regardless of Courant criterion
maxNt = 3; 			%maximum exceedance factor of expected time steps

%perturbation
Bperturb = 2;       %2 branch number
Nperturb = 10;      %10 halfway grid cell number to be perturbed
Zperturb = 0.001;    %0.01 height of perturbation
Tperturb = 100;      %20 time (yr) after which to perturb, allows for cold start

%reports
Nreport  = 10;		%nr of reports for bed and water level curves (>6 cyclic colors)
Trepstep = 10;     %1 default modulo step to store time series data
%location in xcoord for report full time series of base, water, bed levels and Cr
frep     = [ ones(1,length(Topo(:,1))) ]' .*0.05; %fractions of Li for time report

%Tectonics/subsidence in both bifurcates gradual uplift (+) or downwarp (-) of bed
%subsidence 0.03 - 0.05 m / year in literature
etatect = 1.*[0 0 -.05];%0e-3.*[ ones(1,length(Topo(:,1))) ]'; %tectonics rate in m per year
%fraction of L for up/downstream end of 'tectonic block'
ftect   = [ zeros(1,length(Topo(:,1))) ; ones(1,length(Topo(:,1))) ];

%% constants
%Qbreak = Q(1)/25;	%discharge at which one branch is considered closed
If = 1;				%Intermittency
lamp = 0.35;        %Bed Porosity
kc = 1.0; 		    %m Roughness Height, not necessary to calibrate for WhiteColebrook 
bB = 0.5;           %0.6 power on discharge Q for hydraulic geometry
bD = -0.;           %-0.25 power on diameter D for hydraulic geometry
bM = -0.;           %-0.37 power on mud percentage M for hydraulic geometry

%tons/a Imposed annual sediment transport rate fed in from upstream (which must all be carried during floods)
%IF specified then used, ELSE calculated from initial capacity of upstream branch
%Gtf = 1.625E+06;	%

%nodal point constants
k = 1;				%2 POWER OF Wang NODAL POINT RELATION
alb = 1;			%1, r parameter in Ikeda/Bolla NODAL POINT RELATION
%calibrated on Delft3D:
alw = 2;			%3, alpha parameter in Bolla NODAL POINT RELATION
epsilon = 2;		%2, spiral flow calibration parameter
