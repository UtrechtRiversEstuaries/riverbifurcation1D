%preparation and check of surface elevations
%run after initia43

%% check channel profile
Heightspower = repmat(Heights(1,1),1,Nb); 
nS = 2; %power of channel profiles

%titles = {'Biesbosch','Nwe Waterweg','Oude Rijn','Vecht','IJssel'}; %comment away
for telpath = 1:length(seabranches)
    figure(telpath)
    %set(gcf,'units','centimeters','position',[1 1 22 10],'papertype','A4',...
    %    'papertype','A4','paperunits','centimeters','paperposition',[1 1 22 10]);
    %subplot(length(seabranches),1,telpath)
    hold on
    telbranch = seabranches(telpath);
    for telorder = No:-1:2
        plot([xoffset(telbranch) xoffset(telbranch)+L(telbranch)],...
            [Heights(1,telbranch) Heights(2,telbranch)],'.-')
        text(xoffset(telbranch),Heights(1,telbranch)+0.5,...
            Names(telbranch),'color','b')
        text(xoffset(telbranch),-telorder+No,...
            num2str(S(telbranch),'%1.5f'),'color','b')
        %ideal profiles
        Heightspower(telbranch) = min(...
            (Heights(1,1)/(shortest_lengths_to_sea(telpath)^nS)) ...
            * ( shortest_lengths_to_sea(telpath) - xoffset(telbranch) ).^nS , ...
            Heightspower(telbranch)...
            );
        plot(xoffset(telbranch),Heightspower(telbranch),'bo','linewidth',2)
        %if confluence upstream
        if ~isnan(Topo(telbranch,2))
            telbranch2 = Topo(telbranch,2);
            plot([xoffset(telbranch2) xoffset(telbranch2)+L(telbranch2)],...
                [Heights(1,telbranch2) Heights(2,telbranch2)],'r.:')
            text(xoffset(telbranch2),Heights(1,telbranch2)+1,...
                Names(telbranch2),'color','r')
            text(xoffset(telbranch2),-telorder+No,...
                num2str(S(telbranch2),'%1.5f'),'color','r')
            %ideal profiles
            Heightspower(telbranch2) = min(...
                (Heights(1,1)/(shortest_lengths_to_sea(telpath)^nS)) ...
                * ( shortest_lengths_to_sea(telpath) - xoffset(telbranch2) ).^nS , ...
                Heightspower(telbranch2)...
                );
            plot(xoffset(telbranch2),Heightspower(telbranch2),'ro','linewidth',2)
        end
        %continue to next branch
        telbranch = topou{telbranch}(1);
    end
    plot([xoffset(telbranch) xoffset(telbranch)+L(telbranch)],...
        [Heights(1,telbranch) Heights(2,telbranch)],'.-')
    text(xoffset(telbranch),Heights(1,telbranch)+0.5,Names(telbranch))
    text(xoffset(telbranch),-1+No,...
        num2str(S(telbranch),'%1.5f'))
    %ideal profiles
    Heightspower(telbranch) = (Heights(1,1)/(shortest_lengths_to_sea(telpath)^nS)) ...
        * ( shortest_lengths_to_sea(telpath) - xoffset(telbranch) ).^nS;
    plot(xoffset(telbranch),Heightspower(telbranch),'ko','linewidth',2)
    %figure thingies
    plot([0 2e5],[0 0],'k')
    %axis([0 2e5 -0.5 max(Heights(1,:))+0.5]);
    xlabel('distance along river (m)')
    ylabel('height above datum (m)')
    %title([titles(telpath), num2str(shortest_lengths_to_sea(telpath)/1000)]) %comment away
    hold off
end

[[1:Nb]' Heights(1,:)' Heightspower']

%% Check of elevations
% testje
% maxx = 170e3; %k
% testx = [0:maxx/10:maxx]';
% maxz = 13.5;  %m
% power = 2
% multip = maxz/(maxx^power)
% testz = multip.*(maxx-testx).^power;
% plot(testx,testz)

% echte data
% nS = 2;
% aS = Heights(1,1)/(max(shortest_lengths_to_sea)^nS);
% Heightslin = Heights(1,1) - Heights(1,1) .* (xoffset./max(shortest_lengths_to_sea))'; %linear profile
% Heightscur = aS .* ( max(shortest_lengths_to_sea) - xoffset ).^nS;
% plot(xoffset,Heightslin,'.-',xoffset,Heightscur,'.')
% [[1:Nb]' Heights(1,:)' Heightslin' Heightscur']
