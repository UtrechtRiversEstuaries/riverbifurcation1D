%FIGURES
%builds the figures for wb43.m model
%Maarten Kleinhans, February 2008

%%
%report time series
if switchfig1 == 1
    figure
    plot(rep{end}.timeser.timeax(:,1),rep{end}.timeser.baselev(:,1),':')
    hold on
    plot(rep{end}.timeser.timeax(:,1),rep{end}.timeser.waterlev,'-')
    plot(rep{end}.timeser.timeax(:,1),rep{end}.timeser.bedlev,'--')
    title('time series of water levels and bed levels')
    xlabel('time (years)')
    ylabel('height above datum (m)')
    for teltopo = 1:Nb
        text(max(rep{end}.timeser.timeax(:,1)),...
            rep{end}.timeser.waterlev(end,teltopo),num2str(teltopo))
    end
    hold off
end

%%
%report time series
if switchfig2 == 1
    figure
    for telorder = 2:No
        hold on
        plot(rep{end}.timeser.timeax(:,1),...
           rep{end}.timeser.Q(:,Orde(telorder,1):Orde(telorder,2)) ./Sizes(1,1),...
           'linewidth',No+2-telorder)
    end
    title('time series of relative discharge in bifurcated branches')
    xlabel('time (years)')
    ylabel('discharge upstream/downstream (-)')
    axis([0 max(rep{end}.timeser.timeax(:,1)) 0 1]) %Sizes(1,1)
    for teltopo = 2:Nb
        text(max(rep{end}.timeser.timeax(:,1)),rep{end}.timeser.Q(end,teltopo)./Sizes(1,1),num2str(teltopo))
    end
    %legend('Q2/Q1','Q3/Q1','etc',4);
    hold off
end

%%
%report water depths
if switchfig3 == 1
    figure
    for telbifu = 1:nbifurs
        imin = 1+ find( rep{end}.timeser.Q(1,topob{telbifu}(2:3)) == min(rep{end}.timeser.Q(1,topob{telbifu}(2:3))) );
        imax = 1+ find( rep{end}.timeser.Q(1,topob{telbifu}(2:3)) == max(rep{end}.timeser.Q(1,topob{telbifu}(2:3))) );
        if length(imin>1); imin = imin(1); end
        if length(imax>1); imax = imax(end); end
        hold on
        depthmin = rep{end}.timeser.waterlev(:,topob{telbifu}(imin)) - ...
           rep{end}.timeser.bedlev(:,topob{telbifu}(imin));
        depthmax = rep{end}.timeser.waterlev(:,topob{telbifu}(imax)) - ...
           rep{end}.timeser.bedlev(:,topob{telbifu}(imax));
        plot(depthmin,depthmax,'.',depthmin(1),depthmax(1),'o','markersize',5*(nbifurs+1-telbifu))
    end
    title('comparison of water depths in bifurcated branches')
    xlabel('water depth in closing channel (m)')
    ylabel('water depth in growing channel (m)')
end

%%
%report water and bed levels
if switchfig4 == 1
    figure
    hold on
    for teltopo = 1:Nb
       plot(rep{teltopo}.xcor,rep{teltopo}.xii,'--')
       plot(rep{teltopo}.xcor,rep{teltopo}.etai,'-')
    end
    titletxt = [num2str(Durat),'-year evolution of long profile'];
    title(titletxt)
    xlabel('distance along river (m)')
    ylabel('height above datum (m)')
    hold off
end

%%
%report nodal point relation
if switchfig5 == 1
    figure
    for telbifu = 1:nbifurs
        imin = 1+ find( rep{end}.timeser.Q(1,topob{telbifu}(2:3)) == min(rep{end}.timeser.Q(1,topob{telbifu}(2:3))) );
        imax = 1+ find( rep{end}.timeser.Q(1,topob{telbifu}(2:3)) == max(rep{end}.timeser.Q(1,topob{telbifu}(2:3))) );
        if length(imin>1); imin = imin(1); end
        if length(imax>1); imax = imax(end); end
        loglog(rep{end}.timeser.Q(:,topob{telbifu}(imax)) ./ ...
           rep{end}.timeser.Q(:,topob{telbifu}(imin)),...
           rep{end}.timeser.Qsnode(:,topob{telbifu}(imax)) ./ ...
           rep{end}.timeser.Qsnode(:,topob{telbifu}(imin)),'.')
           hold on
    end
    %axis tight
    loglog([.1 10],[.1 10],[.1 10],[.05 20],'--')
    title('Sediment distribution at nodal point')
    xlabel('Q/Q')
    ylabel('Qs/Qs')
    hold off
end

%%
%report width
if (switchfig6 == 1) & (switchBvar == 2)
    figure
    hold on
    for telbifu = 1:nbifurs
       subplot(nbifurs,1,telbifu)
       hold on
       plot(rep{end}.timeser.timeax(:,1),...
          rep{end}.timeser.width1(:,topob{telbifu}(2)), ...
          rep{end}.timeser.timeax(:,1),...
          rep{end}.timeser.width1(:,topob{telbifu}(3)), ...
          rep{end}.timeser.timeax(:,1),...
          (rep{end}.timeser.width1(:,topob{telbifu}(2)) + ...
          rep{end}.timeser.width1(:,topob{telbifu}(3)) ), ...
            'linewidth',2)
          %./ rep{end}.timeser.width1(:,topob{telbifu}(1)), ...
       plot(rep{end}.timeser.timeax(:,1),...
          rep{end}.timeser.width2(:,topob{telbifu}(2)), ...
          rep{end}.timeser.timeax(:,1),...
          rep{end}.timeser.width2(:,topob{telbifu}(3)), ...
          rep{end}.timeser.timeax(:,1),...
          (rep{end}.timeser.width2(:,topob{telbifu}(2)) + ...
          rep{end}.timeser.width2(:,topob{telbifu}(3)) ), ...
          'linewidth',1)
          %./ rep{end}.timeser.width2(:,topob{telbifu}(1)), ...
       ylabel('width2,3,2+3/1 (-)')
        %axis([0 max(rep{end}.timeser.timeax(:,1)) 0 2])
        %legend('B2/B1','B3/B1','(B2+B3)/B1',-1)
       hold off
    end
    xlabel('time (years)')
    title('widths of downstream channels') %Relative 
    hold off
end
