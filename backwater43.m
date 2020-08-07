%backwater solution for use in takke43
Hi{teltopo}(Nx(teltopo)+1) = xi1(teltopo) - etai{teltopo}(Nx(teltopo)+1); %downstream boundary
Ri{teltopo}(Nx(teltopo)+1) = Hi{teltopo}(Nx(teltopo)+1)*Bi{teltopo}(Nx(teltopo)+1)/...
   (2*Hi{teltopo}(Nx(teltopo)+1) + Bi{teltopo}(Nx(teltopo)+1));
for telx = Nx(teltopo):-1:1 %backwater, boundary at Nx+1
   if switchC ==1 %roughness formulation
	   Fbackp = (((etai{teltopo}(telx)-etai{teltopo}(telx+1))/dx) ...
         -( Cf(1) )* ...
         (qf{teltopo}(telx+1)^2/(g* (Hi{teltopo}(telx+1))^3))) ...
	     /(1-(qf{teltopo}(telx+1)^2/(g* (Hi{teltopo}(telx+1))^3)));
	   Hp = Hi{teltopo}(telx+1) - Fbackp*dx;
       Rp = Hp.*Bi{teltopo}(telx)./(2.*Hp + Bi{teltopo}(telx));
	   Fbackc = (((etai{teltopo}(telx)-etai{teltopo}(telx+1))/dx) ...
	     -( Cf(1) )* (qf{teltopo}(telx+1)^2/(g* Hp^3))) ...
	     /(1-(qf{teltopo}(telx+1)^2/(g* Hp^3))); %roughness also adapted!
   elseif switchC ==2
	   Fbackp = (((etai{teltopo}(telx)-etai{teltopo}(telx+1))/dx) ...
         -(( alc2*log10(12.2*Ri{teltopo}(telx+1)/kc) )^(-2))* ...
         (qf{teltopo}(telx+1)^2/(g* (Hi{teltopo}(telx+1))^3))) ...
	     /(1-(qf{teltopo}(telx+1)^2/(g* (Hi{teltopo}(telx+1))^3)));
       Hp = Hi{teltopo}(telx+1) - Fbackp*dx; 
       Rp = Hp.*Bi{teltopo}(telx)./(2.*Hp + Bi{teltopo}(telx));
	   Fbackc = (((etai{teltopo}(telx)-etai{teltopo}(telx+1))/dx) ...
         -(( alc2*log10(12.2*Rp/kc) )^(-2) )* ...
         (qf{teltopo}(telx+1)^2/(g* Hp^3))) ...
	     /(1-(qf{teltopo}(telx+1)^2/(g* Hp^3))); %roughness also adapted!
   end %of switchC
   Hi{teltopo}(telx) = Hi{teltopo}(telx+1) - (1/2)*(Fbackp + Fbackc)*dx;
   Ri{teltopo}(telx) = Hi{teltopo}(telx)*Bi{teltopo}(telx)/...
       (2*Hi{teltopo}(telx) + Bi{teltopo}(telx));
end %backwater
xii{teltopo} = Hi{teltopo}+etai{teltopo};
