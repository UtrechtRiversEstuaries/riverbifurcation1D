%1D morphodynamics model for a network of bifurcated and confluencing rivers
%with Wang and Bolla Pitaluga and new bifurcation nodal point relation
%and variable bifurcate width
%with gradually varied flow and exner
%Maarten Kleinhans, April 2009
%
%uses:
%input44SOMETHING for input parameters (CALL SEPARATELY; REST CALLED BY THIS FILE)
%initia44 for auxiliarly and initial parameters and initiating matrices
%iter43 for flow division at nodal point
%backwater43 for flow
%sedtransp44 for sediment transport
%nodalpt44 for nodal point relations
%transpgrad43 for transport gradients
%bedupd43 for exner bed update including tectonics and bank erosion/deposition
%report44 and 
%endrep44 for output storage
%figs43 for figures (CALL SEPARATELY; REST CALLED BY THIS FILE)

%% initiation
%clear all; close all; clc
%input44SOMETHING
%Auxil. parameters and initialisation grid and bathymetry from down- to upstream:
initia44


%% ****************BEGIN TIME-LOOP*********************************************
tic
waittxt = [num2str(Durat),' yr, ',num2str(max(long)/1000),' km, ',...
    num2str(nbifurs),' bifurs, ',num2str(nconflu),' confl'];
h = waitbar(0,waittxt,'Name','takke44');
telt = 0; 			%time progress
telC = 0;			%keep count of Courant exceedance
timet = 0; 			%actual time
triggerRvar = 0;	%keep track of switching Rvar
triggerSave = 0;	%keep track of saving profiles
timeSave = 0;       %keep track of saving time series
while (timet<Durat & telt<maxNt*Nt ) %& (min(Q)>=Qbreak & telC<10) %DO SOMETHING WITH LOW Q?
    waitbar(timet/Durat,h);
    telt = telt+1;
   
%%	changing boundary conditions per time step
	%Base level rise
	xi1(seabranches) = xi0 + baselevelrise*timet;
	%bend radius (if switchRvar == 1 then constant as specified) 
	if switchRvar == 2 ...
	      & timet-(triggerRvar.*Tbend) >= Tbend
	   triggerRvar = triggerRvar + 1;
	   Rfac = -Rfac;
	elseif switchRvar == 3
        for telR = 1:nbifurs
           x = [timet-dt timet timet+dt] * Lbend(telR) / (2*Tbend(telR));
           y = Abend(telR).*cos(2*pi.*x./Lbend(telR)); %use cos for other phase
           a = sqrt( (x(2)-x(1))^2 + (y(2)-y(1))^2 );
           b = sqrt( (x(3)-x(2))^2 + (y(3)-y(2))^2 );
           c = sqrt( (x(3)-x(1))^2 + (y(3)-y(1))^2 );
           s = (a+b+c)/2;
           warning off
           Radius = sign(y(2))*(a*b*c)/4/sqrt(s*(s-a)*(s-b)*(s-c));
           warning on
           if isnan(Radius) | isempty(Radius)
              Radius = inf;
           end
           Rfac(telR) = Radius/Bi{1}(end);
           %hold on; plot(timet,Rfac,'.')
        end
    end
    %perturbation put in given time after cold start
    if switchperturb == 1
        if timet >= Tperturb
            etai{Bperturb}(Nperturb) = etai{Bperturb}(Nperturb) + Zperturb;
            switchperturb = 0; %switch off so that only perturbed once
        end
    end

%%  BACKWATERS FROM DOWNSTREAM TO UPSTREAM*********************************************
    minimise = repmat(1.1*waterlevelprecision,1,nbifurs); %for iter at bifur
    minimiseold = minimise;
    minimiseoldold = minimise;
    telm = zeros(1,nbifurs); %keep track of number of backwater iterations a bifurcations

    %now iterate between bifurcates in upstream direction and repeat here to conserve discharge
    for telcalliter=1:netiter
        teltopo = 1; telbifur = nbifurs+1; %telconfl = nconflu+1;
        iter43
    end

%     plot(rep{1}.xcor,xii{1},'k',rep{2}.xcor,xii{2},'b',rep{3}.xcor,xii{3},'r--',...
%         rep{4}.xcor,xii{4},'b',rep{5}.xcor,xii{5},'r--',rep{6}.xcor,xii{6},'g:',...
%         rep{7}.xcor,xii{7},'m',rep{7}.xcor,xii{8},'k--')
%     hold on
%     plot(rep{1}.xcor,etai{1},'k',rep{2}.xcor,etai{2},'b',rep{3}.xcor,etai{3},'r--',...
%         rep{4}.xcor,etai{4},'b',rep{5}.xcor,etai{5},'r--',rep{6}.xcor,etai{6},'g:',...
%         rep{7}.xcor,etai{7},'m',rep{8}.xcor,etai{8},'k--')

%%  FROM UPSTREAM TO DOWNSTREAM**************************************************
    %Check activity, compute sediment transport, nodal point relations, gradients and bed updates
    maxCour = 0;
    for teltopo = 1:Nb
       %de-activate active branches below threshold discharge
       if Q(teltopo) < Qbreak & ...
               ( (Activity(1,1) - timet) <= Activity(teltopo,1) & ...
               (Activity(1,1) - timet) >= Activity(teltopo,2) )
           Activity(teltopo,2) = Activity(1,1) - (timet-dt);
       end
       %EXTEND HERE IF NECESSARY TO RE-ACTIVATE

       sedtransp44
       if ( (Activity(1,1) - timet) <= Activity(teltopo,1) )...
               &( (Activity(1,1) - timet) >= Activity(teltopo,2) )
           maxCour = max([maxCour Cour{teltopo}]);
       end
    end
    
    %Compute upstream sediment feed rate at first time step from downstream output
    if ~exist('Gtf','var') & (telt == 1)
        Qsnode(1) = Qi{1}(1); %better estimate of upstreamt transport capacity
%         %or use downstream transport capacity:
%         Qsnodetemp = zeros(size(seabranches));
%         for teltemp = 1:length(seabranches)
%             Qsnodetemp(teltemp) = Qi{seabranches(teltemp)}(1);
%         end
%         Qsnode(1) = sum(Qsnodetemp);
        timeser.Qsnode(telt,1) = Qsnode(1);
    end

    %plot current water surface and bed profiles
    %for teltopo = 1:Nb
    %   hold on; plot(xcoord{teltopo},etai{teltopo},'-'); plot(xcoord{teltopo},xii{teltopo},'r-')
    %end

    %time step adaptation
    dt = dt * CourantCrit / maxCour;
    if dt<dtmin
       dt = dtmin;
    end

    %Nodal point relations
    for telbifur = 1:nbifurs
       nodalpt44
    end
    %Confluences: just sum upstream transport
    for telconfl = 1:nconflu
       Qsnode(topoc{telconfl}(3)) = Qi{topoc{telconfl}(1)}(end) + Qi{topoc{telconfl}(2)}(end);
    end
    %Through-flow nodes: copy upstream transport
    for telthru = 1:nthru
        Qsnode(topot{telthru}(2)) = Qi{topot{telthru}(1)}(end);
    end

    %sediment transport gradients and bed update including tectonics and width-adaptation
    for teltopo = 1:Nb
        if  ( (Activity(1,1) - timet) <= Activity(teltopo,1) )... %then ini in the past
                &( (Activity(1,1) - timet) >= Activity(teltopo,2) ) %and not yet deactivated
            transpgrad43
        else
            dQdx{teltopo} = zeros(Nx(teltopo)+1,1);
        end
        bedupd43 %here to allow for tectonics even if inactive
    end

    %Report
    report44

    timet = timet+dt;
end
%***************END OF TIME-LOOP****************

%% very last storage of results
endrep44

%finish model
close(h);
toc

%figs43