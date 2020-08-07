%nodal point relation for use in takke43

teltopo = topob{telbifur}(1);
teltop2 = topob{telbifur}(2); %FOR POSITIVE Rfac, THIS IS THE outer BEND
teltop3 = topob{telbifur}(3); %FOR POSITIVE Rfac, THIS IS THE inner BEND

if switchnodal == 1 %Wang NODAL POINT relation
   ratioQs23 = ((Bi{teltop2}(1)/Bi{teltop3}(1))^(1-k)) * (Q(teltop2)/Q(teltop3))^k;
   Qsnode(teltop2) = ratioQs23 * (Qi{teltopo}(end)) / (1+ratioQs23); 
   %for convenient comparison to Bolla and Delft3D:
   Qsy(telbifur) = Qsnode(teltop2) - ...
      (Bi{teltop2}(1)/(Bi{teltop2}(1)+Bi{teltop3}(1))) * Qi{teltopo}(end);

elseif switchnodal == 2 %Bolla NODAL POINT relation
   Qy = ( Q(teltop2)-Q(teltop3)-Q(teltopo)*...
      ((Bi{teltop2}(1)-Bi{teltop3}(1))/(Bi{teltop2}(1)+Bi{teltop3}(1))) )/2;
   %H123 = ( ( Hi{teltop2}(1)+Hi{teltop3}(1) )/2 + Hi{teltopo}(end) )/2;
   H123 = Hi{teltopo}(end); %to treat same as mine!
   dzdy = (etai{teltop2}(1) - etai{teltop3}(1))/(Bi{teltopo}(end)/2);
   if switchtransp ~= 3
       qsy = qi{teltopo}(end) * ( (Qy*Hi{teltopo}(end))/(Q(teltopo)*alw*H123) - ...
          dzdy*alb/sqrt(tausi{teltopo}(end)) );
   else
       qsy = qvRsusnode(1) * ...
           qi{teltopo}(end) * ( (Qy*Hi{teltopo}(end))/(Q(teltopo)*alw*H123) - ...
          dzdy*alb/sqrt(tausi{teltopo}(end)) );
   end
   Qsy(telbifur) = qsy*alw*Bi{teltopo}(end);
   Qsnode(teltop2) = Qsy(telbifur) + ...
      (Bi{teltop2}(1)/(Bi{teltop2}(1)+Bi{teltop3}(1))) * Qi{teltopo}(end);

elseif switchnodal == 3 %NEW NODAL POINT relation for meanders
   %Axisymmetric solution for bend profile
   if switchC ==1 %roughness formulation
	   Cf1node = Cf(1); %friction
   elseif switchC ==2
	   Cf1node = ( alc2.*log10(12.2.*Ri{teltopo}(end)./kc) ).^(-2); %friction
   end %of switchC
   if switchtransv == 0
       taus1 = Cf1node.*qf{teltopo}(end)^2./(Rr*g*D(teltopo).*Hi{teltopo}(end).^2); %Shields parameter
       ftheta = 9*(D(teltopo)/Hi{teltopo}(end))^0.3 *taus1^0.5;
   elseif switchtransv == 1
       taus1 = tausi{teltopo}(end);
       ftheta = (1/alb) *taus1^0.5;
   end
   Chezy1 = sqrt( g/Cf1node );
   %DAMPED FOR W DIFFERENCE (means spiral flow reduces to Bolla in the end)
   damppower = 1;
   Bratio = ( min([Bi{teltop2}(1) Bi{teltop3}(1)]) / ...
      max([Bi{teltop2}(1) Bi{teltop3}(1)]) )^damppower;
   A = (2*epsilon/0.4^2) * (1-sqrt(g)/(0.4*Chezy1)) * Bratio; 
   tanbeta = (etai{teltop2}(1) - etai{teltop3}(1))./(Bi{teltopo}(end)/2); %=dzdy
   %sediment distribution
   Qy = ( Q(teltop2)-Q(teltop3)-Q(teltopo)*...
      ((Bi{teltop2}(1)-Bi{teltop3}(1))/(Bi{teltop2}(1)+Bi{teltop3}(1))) )/2;
   V = Qy/(Hi{teltopo}(end)*alw*Bi{teltopo}(end));
   U = Q(teltopo)/(Hi{teltopo}(end)*Bi{teltopo}(end));
   delta = atan(V/U) -atan( A*Hi{teltopo}(end)/(Rfac(telbifur).*Bi{teltopo}(end))); %complete
   tanalpha = ( (sin(delta)-ftheta^-1*tanbeta) /...
      (cos(delta)-ftheta^-1*((etai{teltopo}(end)-etai{teltopo}(end-1))/dx)) ); %approximates 0
   if switchtransp ~= 3
       qsy = qi{teltopo}(end) * tanalpha;
   else
       qsyBED = qvRsusnode(1) * ... does not make much difference in stability
           qi{teltopo}(end) * tanalpha;
       qsySUS = (1-qvRsusnode(1)) * ...
           qi{teltopo}(end) * tan(delta);
       qsy = qsyBED + switchvRsusp * qsySUS;
   end
   Qsy(telbifur) = qsy*alw*Bi{teltopo}(end);
   Qsnode(teltop2) = Qsy(telbifur) + ...
      (Bi{teltop2}(1)/(Bi{teltop2}(1)+Bi{teltop3}(1))) * Qi{teltopo}(end);
   %lambdaR = Hi{teltopo}(end)*Chezy1/sqrt(g); %spiral flow adaptation length (not used yet)

elseif switchnodal == 4 %PARAMETERISED NEW NODAL POINT relation for meanders
   ['not yet implemented']
end; %of nodal point relation

Qsnode(teltop3) = (Qi{teltopo}(end))-Qsnode(teltop2); 
%qsnode2 = Qsnode2/B2(1); qsnode3 = Qsnode3/B3(1);
