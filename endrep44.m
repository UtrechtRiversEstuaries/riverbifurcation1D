%End report for takke44

Nt = telt; Durat = timet;
%profiles
for teltopo = 1:Nb
   rep{teltopo}.xii(:,end)  = xii{teltopo};
   rep{teltopo}.etai(:,end) = etai{teltopo};
	if switchBvar == 2
		rep{teltopo}.width(:,end)	= Bi{teltopo};
	end
end %reporting

%time series
timeSave = timeSave + 1;
timeser.timeax(timeSave,:) 	= [timet,dt];  %time and time step
timeser.baselev(timeSave,:)	= [max(xi0) + baselevelrise*timet,0]; %downstream base level
timeser.Q(timeSave,:)		= Q; %discharge
timeser.Qsnode(timeSave,:)	= Qsnode; %sediment input
timeser.Qsy(timeSave,:)		= Qsy; 	%transverse sediment flux at bifur
for teltopo = 1:Nb
   timeser.waterlev(timeSave,teltopo)  = xii{teltopo}(irep{teltopo});
   timeser.bedlev(timeSave,teltopo)    = etai{teltopo}(irep{teltopo});
   if switchBvar == 2                    %when width is varying
      timeser.width1(timeSave,teltopo) = Bi{teltopo}(1);
      timeser.width2(timeSave,teltopo) = Bi{teltopo}(end);
      timeser.Twidth(timeSave,teltopo) = Tw{teltopo}(1);
   end
   timeser.shields(timeSave,teltopo)   = tausi{teltopo}(irep{teltopo});
end

%remove redundant time series length
timeser.timeax(any(isnan(timeser.timeax)'),:) = [];
timeser.baselev(any(isnan(timeser.baselev)'),:) = [];
timeser.waterlev(any(isnan(timeser.waterlev)'),:) = [];
timeser.bedlev(any(isnan(timeser.bedlev)'),:) = [];
timeser.Q(any(isnan(timeser.Q)'),:) = [];
timeser.Qsnode(any(isnan(timeser.Qsnode)'),:) = [];
if nbifurs == 1
    timeser.Qsy(isnan(timeser.Qsy)) = [];
else
    timeser.Qsy(any(isnan(timeser.Qsy)'),:) = [];
end
if switchBvar == 2 %when width is varying
	timeser.width1(any(isnan(timeser.width1)'),:) = [];
	timeser.width2(any(isnan(timeser.width2)'),:) = [];
	timeser.Twidth(any(isnan(timeser.Twidth)'),:) = [];
end
timeser.shields(any(isnan(timeser.shields)'),:) = [];

%add timeseries structure to rep cell array as number_of_branches+1^th cell
rep{Nb+1}.timeser = timeser;
