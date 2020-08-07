%Report for takke44

%time series
if telt == 1
    timeser.timeax(telt,:) 	= [timet,dt];  %time and time step
    timeser.baselev(telt,:)	= [max(xi0) + baselevelrise*timet,0];  	%downstream base level and possibly other param
    timeser.Q(telt,:)		= Q; %discharge
    timeser.Qsnode(telt,:)	= Qsnode; %sediment input
    timeser.Qsy(telt,:)		= Qsy; 	%transverse sediment flux at bifur
    for teltopo = 1:Nb
       timeser.waterlev(telt,teltopo) = xii{teltopo}(irep{teltopo});
       timeser.bedlev(telt,teltopo)   = etai{teltopo}(irep{teltopo});
       if switchBvar == 2 %when width is varying
          timeser.width1(telt,teltopo) = Bi{teltopo}(1);
          timeser.width2(telt,teltopo) = Bi{teltopo}(end);
          timeser.Twidth(telt,teltopo) = Tw{teltopo}(1);
       end
       timeser.shields(telt,teltopo)   = tausi{teltopo}(irep{teltopo});
    end
elseif floor(telt/Trepstep) >= timeSave
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
end %reporting time series    

%profiles
if telt == 1
   for teltopo = 1:Nb
      rep{teltopo}.xii(:,1)  = xii{teltopo};
      rep{teltopo}.etai(:,1) = etai{teltopo};
      if switchBvar == 2
          rep{teltopo}.width(:,1)	= Bi{teltopo};
      end
   end
elseif timet-(triggerSave*Durat/Nreport) >= Durat/Nreport
   triggerSave = triggerSave + 1;
   for teltopo = 1:Nb
      rep{teltopo}.xii(:,triggerSave + 1)  = xii{teltopo};
      rep{teltopo}.etai(:,triggerSave + 1) = etai{teltopo};
	  if switchBvar == 2
		  rep{teltopo}.width(:,triggerSave + 1)	= Bi{teltopo};
	  end
   end
end %reporting
