%iterate backwaters through the network for use in takke43

telo = No;
while telo>0 %work from downstream orders to upstream
    teltopo = Orde(telo,1);
    telbifurold = telbifur;
    while teltopo<=Orde(telo,2) %work from low to high branch numbers

        %most upstream branch
        if isnan( Topo(teltopo,1) )
            %if true then no upstream channel exists at all so this is the most upstream channel
            %do one backwater, discharge remains constant
            backwater43 %for this teltopo


        %bifurcation
        elseif ~isnan( Topo( Topo(teltopo,1), 4) ) 
           %if true then there exists a second channel attached to the upstream channel of the present channel
           %do iterate backwaters in bifurcation to obtain discharge in both branches
           telbifur = telbifur - 1;
           teltopoold = teltopo;
           teltopob = [teltopo Topo( Topo(teltopo,1), 4)];
           tempteltop2 = teltopob(find( [Bi{teltopob(1)}(1) Bi{teltopob(2)}(1)] ...
               == min([Bi{teltopob(1)}(1) Bi{teltopob(2)}(1)]) ));
           if length(tempteltop2)==1; %in case equal branches
               teltop2(telbifur)=tempteltop2;
               teltop3(telbifur)=teltopob(find(teltopob~=tempteltop2));
           else
               teltop2(telbifur)=tempteltop2(1);
               teltop3(telbifur)=teltopob(find(teltopob~=tempteltop2(1)));
           end

           %update discharge with upstream changes
           Q(teltop2(telbifur)) = Q(teltop2(telbifur)) + ...
               ( Bi{teltop2(telbifur)}(1) / (Bi{teltop2(telbifur)}(1)+Bi{teltop3(telbifur)}(1)) ) .* ...
               ( Q(Topo(teltopo,1)) - Qold(Topo(teltopo,1)) );
           Q(teltop3(telbifur)) = Q(Topo(teltopo,1))-Q(teltop2(telbifur)); %mass conservation
           qf{teltop2(telbifur)} = Q(teltop2(telbifur))./Bi{teltop2(telbifur)}; %m^2/s specific discharge
           qf{teltop3(telbifur)} = Q(teltop3(telbifur))./Bi{teltop3(telbifur)};
           %remember discharge
           Qold(teltop2(telbifur))=Q(teltop2(telbifur)); 
           Qoldold(teltop2(telbifur))=Qold(teltop2(telbifur));
           Qold(teltop3(telbifur))=Q(teltop3(telbifur)); 
           Qoldold(teltop3(telbifur))=Qold(teltop3(telbifur));
           
           %now loop to get same water level and the concurrent discharge division
           minimiseold(telbifur)=minimise(telbifur); 
           minimiseoldold(telbifur)=minimiseold(telbifur);
           while ( (abs(minimise(telbifur))>waterlevelprecision) ...
                 & (telm(telbifur)<maxiter) )
              telm(telbifur) = telm(telbifur)+1;
              teltopo = teltop2(telbifur);
              backwater43
              teltopo = teltop3(telbifur);
              backwater43
              %now continue with equating water heights in branches 2 and 3
              minimise(telbifur) = xii{teltop2(telbifur)}(1)-xii{teltop3(telbifur)}(1);
              %THE NEW GUESS FOR DISCHARGE:
              Q2temp{telbifur} = Q(teltop2(telbifur))* ...
                 ((Hi{teltop2(telbifur)}(1)- 1*minimise(telbifur) )/Hi{teltop2(telbifur)}(1));
              if telm(telbifur)>2 & sign(minimiseold(telbifur))~=sign(minimise(telbifur))
%               if telm(telbifur)>2 & ...
%                       ( sign(minimiseold(telbifur))~=sign(minimise(telbifur)) | ...
%                       sign(minimiseoldold(telbifur))~=sign(minimiseold(telbifur)) )
                   Q(teltop2(telbifur)) = (Qold(teltop2(telbifur)) + Q2temp{telbifur} )/2;
%                    Q(teltop2(telbifur)) = ...
%                        (Qoldold(teltop2(telbifur)) + Qold(teltop2(telbifur)) + Q2temp{telbifur} )/3;
              else
                   Q(teltop2(telbifur)) = Q2temp{telbifur};
              end
              Q(teltop3(telbifur)) = Q(Topo(teltopo,1))-Q(teltop2(telbifur)); %mass conservation
              qf{teltop2(telbifur)} = Q(teltop2(telbifur))./Bi{teltop2(telbifur)}; %m^2/s specific discharge
              qf{teltop3(telbifur)} = Q(teltop3(telbifur))./Bi{teltop3(telbifur)};
           end
           xi1(Topo(teltopo,1)) = ( xii{teltop2(telbifur)}(1) + xii{teltop3(telbifur)}(1) )/2;
           teltopo = teltopoold + 1;


        %through-flow (for getting order of upstream branches for confluence right)
        elseif isnan( Topo(teltopo,2) )
            %if true then there exists neither a second upstream channel (confluence) nor a bifurcate (previous if)
            Qoldold(teltopo) = Qold(teltopo); 
            Qold(teltopo) = Q(teltopo); 
            Q(teltopo) = Q(Topo(teltopo,1)); %mass conservation
            qf{teltopo} = Q(teltopo)./Bi{teltopo}; %m^2/s specific discharge
            %do one backwater
            backwater43 %for this teltopo
            %make upstr water lev the downstr bound for the upstr channel
            xi1(Topo(teltopo,1)) = xii{teltopo}(1);


        %confluence
        elseif ~isnan( Topo(teltopo,2) )
            %if true then a second upstream channel exists (confluence)
            Qoldold(teltopo) = Qold(teltopo); 
            Qold(teltopo) = Q(teltopo); 
            Q(teltopo) = Q(Topo(teltopo,1)) + Q(Topo(teltopo,2)); %mass conservation
            qf{teltopo} = Q(teltopo)./Bi{teltopo}; %m^2/s specific discharge
            %do backwater with same downstream water level as second upstream channel
            backwater43 %for this teltopo
            %make upstr water lev the downstr bound for two upstr channels
            xi1(Topo(teltopo,1)) = xii{teltopo}(1);
            xi1(Topo(teltopo,2)) = xii{teltopo}(1);

        end

        teltopo = teltopo+1;

    end %of index for branch

    telo = telo-1;

end %of index for branch order
