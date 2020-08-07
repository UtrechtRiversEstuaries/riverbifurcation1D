%calculates transport rates for takke44
qsi{teltopo} = zeros(Nx(teltopo)+1); 

if switchtransp == 1 %choice transport equation EH
   if switchC ==1 %roughness formulation
	   Cfi{teltopo} = Cf1; %friction
   elseif switchC ==2
	   Cfi{teltopo} = ( alc2.*log10(12.2.*Ri{teltopo}./kc) ).^(-2); %friction
   end %of switchC
   tausi{teltopo} = Cfi{teltopo}.*qf{teltopo}.^2./(Rr*g*D(teltopo).*Hi{teltopo}.^2); %Shields parameter
   qsi{teltopo} = (alt./Cfi{teltopo}).*tausi{teltopo}.^nt; %EH, Einstein par

elseif switchtransp == 2 %MPM
   if switchCp == 1
       Cfi{teltopo} = ( alc2.*log10(12.2.*Ri{teltopo}./kc) ).^(-2); %total friction
   else
       Cfi{teltopo} = ( alc2.*log10(12.2.*Ri{teltopo}./(2.5*D(teltopo))) ).^(-2); %grain friction
   end
   tausi{teltopo} = Cfi{teltopo}.*qf{teltopo}.^2./(Rr*g*D(teltopo).*Hi{teltopo}.^2); %Shields parameter
   qstemp = zeros(size(tausi{teltopo}));
   transp = find(tausi{teltopo}>tausc); %only for above motion MPM
   qstemp = alt.*(tausi{teltopo}(transp)-tausc).^nt; %MPM, Einstein par
   qsi{teltopo} = qstemp; %MPM, Einstein par
   if length(transp)<length(qstemp)
       ['attention! Sediment below motion in branch ' num2str(teltopo) ...
           ' at ' num2str(timet) ' year (' num2str(Activity(1,1)-timet) ' BP)']
   end

elseif switchtransp == 3 %van Rijn BEDLOAD AND SUSPENDED LOAD
   if switchCp == 1
       Cfi{teltopo} = ( alc2.*log10(12.2.*Ri{teltopo}./kc) ).^(-2); %total friction
       za = kc;
   else
       Cfi{teltopo} = ( alc2.*log10(12.2.*Ri{teltopo}./(2.5*D(teltopo))) ).^(-2); %grain friction
       za = 2.5*D(1);
   end
   tausi{teltopo} = Cfi{teltopo}.*qf{teltopo}.^2./(Rr*g*D(teltopo).*Hi{teltopo}.^2); %Shields parameter
   T = (tausi{teltopo} - tausc)./tausc;
   Dstar = D(teltopo)*(Rr*g/1.2e-6^2).^(1/3); %Bonnefille dimensionless grain size
   qstempBEDL = zeros(size(tausi{teltopo})); 
   transp = find(tausi{teltopo}>tausc); %only for above motion
   if length(transp)<length(qstempBEDL)
       ['attention! Sediment below motion in branch ' num2str(teltopo) ...
           ' at ' num2str(timet) ' year (' num2str(Activity(1,1)-timet) ' BP)']
   end
   qstempBEDL(transp) = alt.*T(transp).^nt.*Dstar^-0.3; %van Rijn BEDLOAD, Einstein par
   qstempSUSP = zeros(size(qstempBEDL)); ca = qstempSUSP; %SUSPENDED
   veloc = Q(teltopo)./Hi{teltopo}./Bi{teltopo};
   ca(transp) = als.*(D(teltopo)/za).*T(transp).^nt.*Dstar^-0.3;
   ca(ca>(1-lamp)) = (1-lamp);
   if switchC ==1 %roughness formulation, TOTAL shear stress for u* needed
	   Cfi{teltopo} = Cf1; %friction
   elseif switchC ==2
	   Cfi{teltopo} = ( alc2.*log10(12.2.*Ri{teltopo}./kc) ).^(-2); %friction
   end %of switchC
   tausi{teltopo} = Cfi{teltopo}.*qf{teltopo}.^2./(Rr*g*D(teltopo).*Hi{teltopo}.^2); %Shields parameter
   ustar = sqrt(tausi{teltopo}.*Rr*g*D(teltopo));
   %ws1 = (1.2e-6./D(teltopo)).*(sqrt(10.36^2+1.049*(1-0)^4.7.*Dstar.^3)-10.36); 
   ws1 = 1.1.*sqrt(Rr*g*D(teltopo));
   beta = 1+2.*(ws1./ustar).^2; 
   beta(beta>2) = 2;
   Z = ws1./(beta.*0.4.*ustar);
   F = ((za./Hi{teltopo}).^Z-(za./Hi{teltopo}).^1.2)./(((1-za./Hi{teltopo}).^Z).*(1.2-Z));
   qstempSUSP = (F.*veloc.*za./Hi{teltopo}.*ca)./(sqrt(Rr*g*D(teltopo))*D(teltopo)); %van Rijn SUSPENDED
   qsi{teltopo} = qstempBEDL + switchvRsusp * qstempSUSP; %van Rijn TOTAL, Einstein par
   qvRsusnode(teltopo) = qstempBEDL(end) / (qstempBEDL(end) ...
       + switchvRsusp * qstempSUSP(end));
   
end; %choice transport equation
qi{teltopo} = qsi{teltopo}.*sqrt(Rr*g*D(teltopo))*D(teltopo); %m^2/s
Qi{teltopo} = qi{teltopo}.*Bi{teltopo}; %m^3/s

%Courant numbers for bed disturbance celerities
cbed{teltopo} = max(2*nt.* (qi{teltopo}./Hi{teltopo}./(1-lamp)) );
Cour{teltopo} = cbed{teltopo}*dt*annual/dx;
if ( (Activity(1,1) - timet) <= Activity(teltopo,1) )...
        &( (Activity(1,1) - timet) >= Activity(teltopo,2) )
    if Cour{teltopo}>1
       telC = telC+1;
       ['attention! Courant=' num2str(Cour{teltopo}) ' in branch ' num2str(teltopo) ...
           ' at ' num2str(timet) ' year (' num2str(Activity(1,1)-timet) ' BP)']
    end
end