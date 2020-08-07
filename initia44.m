% Prepares matrices for takke44

% Auxiliarly parameters, constants and initiation of matrices to start the model takke43
Nt = Durat/dt;		%Number of time steps (approximate for variable time stepping)

%% extra switches mostly for perturbations module
if exist('switchCp','var') ~= 1
    switchCp = 2;
end

if exist('switchvRsusp','var') ~= 1
    switchvRsusp = 1;
end

if exist('switchtransv','var') ~= 1
    switchtransv = 1;
end

if exist('switchperturb','var') ~= 1
    switchperturb = 0;
end

%% constants for sediment transport
g = 9.81;			%9.81 earth, 3.74 Mars m/s^2 gravitational acceleration
alc1 = (sqrt(g)/log10(exp(1))/0.4); %18.02 roughness constant for Colebrook-white
alc2 = (1/log10(exp(1))/0.4); %5.76 roughness constant for Colebrook-white
Rr = 1.65;			%submerged specific gravity of sediment
annual = 60*60*24*365.25; %number of seconds in a year

if switchtransp == 1 %choice transport equation
	alt = 0.05;		%coefficient in Engelund-Hansen total bed material type relation
	nt = 2.5;		%2.5exponent in Engelund-Hansen total bed material type relation
elseif switchtransp == 2
    dt = dt/3;
    CourantCrit = CourantCrit/3;
    dtmin = dtmin/3;
    maxNt = 3*maxNt;
	alt = 8;		%5.7 coefficient in Meyer-Peter Muller type load relation
    nt = 1.5;		%1.5 exponent in Meyer-Peter Muller type load relation
    tausc = 0.047;	%0.03 critical Shields stress in Meyer-Peter Muller type load relation
elseif switchtransp == 3
    dt = dt/3;
    CourantCrit = CourantCrit/3;
    dtmin = dtmin/3;
    maxNt = 3*maxNt;
    als = 0.015;    %0.015 coefficient in van Rijn reference concentration suspended load
    alt = 0.1;		%0.1 coefficient in van Rijn bed load relation
    nt = 1.5;		%1.5 exponent in van Rijn bed load relation
    tausc = 0.047;	%0.047 critical Shields stress in van Rijn bed load relation
    
end %choice transport equation

%% Network topology
%get old matrices temporarily
Q = Sizes(:,1)'; %m^3/s Flood discharge
L = Sizes(:,2)'; %length of the branches
D = Sizes(:,3)'./1000; %grain sizes per branch in m
M = Sizes(:,4)'; %mud percentage per branch

%matrix sizes
Nb = length(Topo(:,1));
No = length(Orde(:,1));
%number of iterations in flow network in case not specified
if exist('netiter','var') ~= 1
    netiter = No;
end

%determine which branches have free downstream nodes (at sea)
seabranches = find(isnan(Topo(:,3))==1);

%make old topou and topod matrices
%upstream connections: (none if upstream boundary, one if bifur, two if confluence)
%downstream connections: (none if sea, one if confluence, two if bifurcation)
for tel=1:Nb
    topou{tel}=[ Topo(tel,find(~isnan(Topo(tel,1:2)))) ];
    topod{tel}=[ Topo(tel,2+find(~isnan(Topo(tel,3:4)))) ];
end
%topou{1}=[]; topou{2}=[1]; topou{3}=[1]; topou{4}=[2]; topou{5}=[2]; topou{6}=[3]; topou{7}=[3];
%topod{1}=[2 3]; topod{2}=[4 5]; topod{3}=[6 7]; topod{4}=[]; topod{5}=[]; topod{6}=[]; topod{7}=[];

%define downstream conditions for backwater; number branches; empty for sea level
%  e.g. 4(6)->31(4)+32(5) for confluence, 31+32->3(3) for bifur, 
%  2(2)+3->1(1) for nodal bifur:
%bifurcations: up and downstream branch numbers
%confluences: up and downstream branch numbers

bifuri = find(~isnan(Topo(:,4))); %branches that have bifurcation DOWNSTREAM
nbifurs = length(bifuri);	%number of bifurcations
for tel=1:nbifurs
    %topob{length(bifuri)-tel+1}=[ bifuri(tel) Topo(bifuri(tel),3:4) ];
    topob{tel}=[ bifuri(tel) Topo(bifuri(tel),3:4) ];
end
confli = find(~isnan(Topo(:,2))); %branches that have confluence UPSTREAM
nconflu = length(confli);	%number of confluences
if nconflu>0
    for tel=1:nconflu
        topoc{tel}=[ Topo(confli(tel),1:2) confli(tel) ];
    end
else
    topoc=[];
end

%branches that have throughflow node UPSTREAM
nthru = length(Connections);	%number of through-flow nodes (simple connections to get order right)
%nthru = length(thrui);	%number of through-flow nodes (simple connections to get order right)
if nthru>0
    upstrdone = [1];
    downstrdone = [seabranches'];
    for tel=1:nbifurs
        upstrdone = [upstrdone topob{tel}(2:3)];
    end
    for tel=1:nconflu
        downstrdone = [downstrdone topoc{tel}(1:2)];
    end
    %now find whether done as other branch
    for tel=1:length(upstrdone)
        downstrdone = [downstrdone topou{upstrdone(tel)}];
    end
    for tel=1:length(downstrdone)
        upstrdone = [upstrdone topod{downstrdone(tel)}];
    end
    upstrnotdone = setdiff([1:Nb] , upstrdone);
    downstrnotdone = setdiff([1:Nb] , downstrdone);

    %thrui = setdiff( find( isnan(Topo(:,2)) & isnan(Topo(:,4)) ) , ...
    %   [Topo(bifuri,3); Topo(bifuri,4); Topo(confli,2); Topo(confli,2) ] );
    for tel=1:nthru
%        topot{tel}=[ Topo(thrui(tel),1) thrui(tel) ];
        topot{tel}=[downstrnotdone(tel) upstrnotdone(tel) ];
    end
else
    topot=[];
end


%% Bends
%Upstream bend: relative bend radius R/W at each bifurcation
Rfac = Bifurcations(2:end,2)';
%bifurcate 2 is outer bend for positive (for variable R switchRvar=2)
Tbend = Bifurcations(2:end,3)';
%period (yr) at which bend radius changes sign, REAL period is 2Tbend!
							%(for variable R switchRvar=2 and 3 resp)
Lbend = Bifurcations(2:end,4)';
%wave length of sinusoidal double bend (for switchRvar=3)
Abend = Bifurcations(2:end,5)';
%bend amplitude (for switchRvar=3)

% Nb = 1 + 2*nbifurs +nconflu; %number of branches
teltop2 = repmat(NaN,1,nbifurs);
teltop3 = teltop2;
Nx = L./dx;		%number of spatial steps (=nodes-1) (excluding ghost node)
if any(Nx-round(Nx)~=0)
   ['no integer number of spatial steps']
end


%% Elevations
%sea level
xi1 = repmat(NaN,1,Nb); %downstream boundary
xi1(seabranches) = xi0; %impose sea level

%upstream and downstream bed levels for each branch
Heights = repmat(NaN,2,Nb);
%Heights(1,1) = Height(1);
Heights(1,1) = Bifurcations(1);
Heights(2,seabranches) = xi0;
for teltopo=1:nbifurs %bifurcations
   Heights(2,topob{teltopo}(1))   = Bifurcations(teltopo+1,1);
   Heights(1,topob{teltopo}(2:3)) = Bifurcations(teltopo+1,1);
end
for teltopo=1:nconflu %confluences
   Heights(2,topoc{teltopo}(1:2)) = Confluences(teltopo);
   Heights(1,topoc{teltopo}(3))   = Confluences(teltopo);
end
for teltopo=1:nthru
   Heights(2,topot{teltopo}(1)) = Connections(teltopo);
   Heights(1,topot{teltopo}(2)) = Connections(teltopo);
end

%channel gradients:
S = ( Heights(1,:) - Heights(2,:) )./L;

%test topology: how long are different pathways to the sea
%cumulative length with each branch
long = zeros(Nb,1);
for teltopo=1:Nb
   if length(topou{teltopo}) == 1       %bifurcation
      long(teltopo) = long(topou{teltopo}) + L(teltopo);
   elseif length(topou{teltopo}) == 2   %confluence
      long(teltopo) = max(long(topou{teltopo})) + L(teltopo);
   elseif length(topou{teltopo}) == 0   %upstream boundary
      long(teltopo) = L(teltopo);
   end
end
xoffset = long-L'; %for reported x-coordinates
shortest_lengths_to_sea = long(seabranches);

%% Estimation of ambient river conditions WITH NORMAL FLOW APPROXIMATION
aB = Bup/(Q(1)^bB * D(1)^bD * M(1)^bM);  %constant for hydraulic geometry width predictor
B = aB .*Q.^bB .*D.^bD .*M.^bM;        %m Channel Width, 550 for Rijn 1792

%just initial guesses
u = repmat(1,1,Nb);
H = (Q.^2 .*kc^(1/3) ./(8.1^2)./g./ S./B.^2).^(3/10); %convenient Parker equation
Htemp = 0.9.*H;
R = H.*B./(2.*H+B);
%teltemp = 0;
while abs( max(H-Htemp) ) > waterlevelprecision
   C1 = alc1.*log10((12.2.*R)./kc);
   if switchC ==1 %roughness formulation
      C1 = repmat(C1(1),size(C1));
   end %of switchC
   u = C1.*sqrt(R.*S);
   Htemp = Q./(u.*B);
   H = ( 2.*H+Htemp )./3;
   R = H.*B./(2.*H+B);
   %hold on; plot(teltemp,H,'.',teltemp,Htemp,'o');
   %teltemp = teltemp+1;
end
Cf = g./( C1.^2 );
u = Q./B./H;
Qini = Q;
Qold = Q;
Qoldold = Q;
if ~exist('Qbreak','var')
    Qbreak = Sizes(1,1)/25;
end
%Hmin = min(H)/2;

%initial channel depths
eta0 = max(xi0) - H;
%Frtest = (Q./B./H)./sqrt(g.*H);

%% Estimation of transport rate in UPSTREAM channel
if switchtransp == 1 %EH; choice transport equation
   Cf1 = Cf(1);
%   Cf1a = ( 5.76.*log10(12.2.*R(1)./kc) ).^(-2); %friction
   taus1 = Cf1.*Q(1)^2./(Rr*g*D(1).*H(1).^2.*B(1).^2); %Shields parameter
   qs1 = (alt./Cf1).*taus1.^nt; %EH, Einstein par

elseif switchtransp == 2 %MPM
   if switchCp == 1
       Cf1 = ( alc2.*log10(12.2.*R(1)./kc) ).^(-2); %total friction
   else
       Cf1 = ( alc2.*log10(12.2.*R(1)./(2.5*D(1))) ).^(-2); %grain friction
   end
   taus1 = Cf1.*Q(1)^2./(Rr*g*D(1).*H(1).^2.*B(1).^2); %Shields parameter
   if taus1>tausc %only for above motion MPM
       qs1 = alt.*(taus1-tausc).^nt; %MPM, Einstein par
   else
       qs1 = 0;
       ['attention! Sediment below motion at model start']
   end %above motion

elseif switchtransp == 3 %van Rijn
   if switchCp == 1
       Cf1 = ( alc2.*log10(12.2.*R(1)./kc) ).^(-2); %total friction
       za = kc;
   else
       Cf1 = ( alc2.*log10(12.2.*R(1)./(2.5*D(1))) ).^(-2); %grain friction
       za = 2.5*D(1);
   end
   taus1 = Cf1.*Q(1)^2./(Rr*g*D(1).*H(1).^2.*B(1).^2); %Shields parameter
   T = (taus1 - tausc)./tausc;
   ustar = sqrt(taus1*Rr*g*D(1));
   Dstar = D(1)*(Rr*g/1.2e-6^2).^(1/3); %Bonnefille dimensionless grain size
   ws1 = (1.2e-6/D(1)).*(sqrt(10.36^2+1.049*(1-0)^4.7.*Dstar.^3)-10.36); 
   if taus1>tausc %only for above motion MPM
       qs1bed = alt.*T.^nt.*Dstar^-0.3; %van Rijn, Einstein par
       ca = als*(D(1)/za)*T^nt*Dstar^-0.3;
       ca(ca>(1-lamp)) = (1-lamp);
       beta = 1+2.*(ws1./ustar).^2; 
       beta(beta>2) = 2;
       Z = ws1./(beta.*0.4.*ustar);
       F = ((za./H(1)).^Z-(za./H(1)).^1.2)./(((1-za./H(1)).^Z).*(1.2-Z));
       qs1sus = F.*u(1).*za./H(1).*ca./sqrt(Rr*g*D(1))*D(1); %van Rijn SUSPENDED
       qs1 = qs1bed + switchvRsusp * qs1sus;
   else
       qs1 = 0;
       ['attention! Sediment below motion at model start']
   end %above motion

end; %choice transport equation
qt	= qs1*sqrt(Rr*g*D(1))*D(1); %m^2/s Volume sediment transport rate per unit width (at flood)
Gt	= qt*B(1)*annual*(Rr+1)*If; %tons/a Ambient annual sediment transport rate in tons per annum (averaged over entire year)
%Calculation of ultimate conditions imposed by a modified rate of sediment input
%m^2/s Upstream imposed volume sediment transport rate per unit width (at flood)
%AND NOW Specification of input sediment transport m^2/s (during floods) at GHOST NODE
qsnode = repmat(NaN,1,Nb);
if exist('Gtf','var')
    qtG = Gtf/(Rr+1)/annual/B(1)/If;
    qsnode(1) = qtG; %qtG=feed
    Qsnode = B.*qsnode;
else
    qsnode(1) = qt;  %qt=capacity
	Qsnode = B.*qsnode;
    %will be replaced by better calculation in takke43 based on backwater flow
end

%% initialisation of matrices
Qsy = repmat(NaN,1,nbifurs); %cross transport at bifurcation
if switchtransp == 3
    qvRsusnode = zeros(Nb,1); %make matrix for suspended load for bifurcations
end

for teltopo = 1:Nb
   xcoord{teltopo} = (0:dx:L(teltopo))'; %make base grid
   irep{teltopo} = round(frep(teltopo)*Nx(teltopo)); %position for which output is reported
   if irep{teltopo} == 0
      irep{teltopo} = 1;
   end
   itect{teltopo} = find( (xcoord{teltopo}>round(ftect(1,teltopo)*L(teltopo))) &...
      (xcoord{teltopo}<round(ftect(2,teltopo)*L(teltopo))) ); %position at which tectonics is applied
   Bi{teltopo} = repmat(B(teltopo),length(xcoord{teltopo}),1); %width
   qf{teltopo} = Q(teltopo)./Bi{teltopo}; %m^2/s specific discharge
   dQdx{teltopo} = zeros(size(xcoord{teltopo}));
   xii{teltopo} = repmat(NaN,length(xcoord{teltopo}),1); %local water level
   etai{teltopo} = Heights(1,teltopo) -S(teltopo).*xcoord{teltopo} + eta0(teltopo); %incl correct inlet jump
   %OUD: etai = S.*(max(xcoord)-xcoord) + xii0 - H(teltopo); %incl correct inlet jump
   %storage matrix
   rep{teltopo}.xcor = xcoord{teltopo} + xoffset(teltopo);
   rep{teltopo}.xii  = repmat(NaN,length(xcoord{teltopo}),Nreport+1);
   rep{teltopo}.etai = rep{teltopo}.xii;
   if switchBvar == 2
       rep{teltopo}.width = rep{teltopo}.xii;
   end
   %hold on; plot(xcoord{teltopo},etai{teltopo},'-')
end
Hi = xii; Ri = xii; Tw = xii; Bprevious = Bi; etaitemp = etai;

% figure; hold on
% for teltopo = 1:Nb
%    plot(xcoord{teltopo},etai{teltopo},'-')
%    pause
% end

%make storage matrix for time series; later add to rep cell array
if exist('Trepstep','var') ~= 1
    Trepstep = 1;
end
timeserlength       = round(maxNt*Nt/Trepstep);
timeser.timeax      = repmat(NaN,timeserlength,2);  %time and time step
timeser.baselev     = repmat(NaN,timeserlength,2);  %downstream base level and possibly other param
timeser.waterlev	= repmat(NaN,timeserlength,Nb); %water level at report location
timeser.bedlev		= repmat(NaN,timeserlength,Nb); %bed level at report location
timeser.Q			= repmat(NaN,timeserlength,Nb); %discharge
timeser.Qsnode		= repmat(NaN,timeserlength,Nb); %sediment input
timeser.Qsy			= repmat(NaN,timeserlength,nbifurs); %transverse sediment flux at bifur
if switchBvar == 2    %when width is varying
   timeser.width1	= repmat(NaN,timeserlength,Nb); %upstream width
   timeser.width2	= repmat(NaN,timeserlength,Nb); %downstream width
   timeser.Twidth	= repmat(NaN,timeserlength,Nb); %time scale for width (upstream)
end
timeser.shields		= repmat(NaN,timeserlength,Nb); %sediment input

%preparation of Activity matrix if not specified (older versions)
if exist('Activity','var') == 0
    Activity = repmat([Durat 0],Nb,1);
end


%% Unnecessary
%Estimation ultimate equilibrium Shields number, NOT USED but same as in Parker spreadsheet
%tausu = (tausc+(qtG/sqrt(Rr*g*D(1))/D(1)/alt)^(1/nt));
%Ultimate slope to which the bed must aggrade
%Su = ((Rr*D(1)*tausu)^(10/7))*((8.1^2*B(1)^2*g/Q(1)^2/kc^(1/3)))^(3/7);
%m Ultimate flow depth (at flood)
%Hu =Rr*D(1)*tausu/Su;

%Specification of Imposed Downstream Water Surface Elevation, NOT USED but same as in Parker spreadsheet
%The user imposes a water surface elevation xi0.
%Frni = sqrt(Q(1)^2/g/B(1)^2/H(1)^3); %Initial normal Froude number; must be < 1 to proceed
%Frnu = sqrt(Q(1)^2/g/B(1)^2/Hu^3); %Ultimate normal Froude number; must be < 1 to proceed
%Hc1 = (Q(1)^2/B(1)^2/g)^(1/3); %m Critical depth (FROUDE=1)
%ximin = Hc1-eta0; %m Minimum possible downstream water surface elevation; compare to xim
