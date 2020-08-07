%bed update including tectonics and width-adaptation for use in takke43

if switchBvar == 1
	etai{teltopo} = etai{teltopo} - (1/(1-lamp)).*(dQdx{teltopo}./Bi{teltopo}).*If*annual*dt; %Exner
	etai{teltopo}(itect{teltopo}) = etai{teltopo}(itect{teltopo}) + etatect(teltopo)*dt;
elseif switchBvar == 2 %bank erosion/deposition, adaptation effects included!
%	Tw = 1;  								%1. CONSTANT ADAPTATION
%	Tfac = 1;								%2. Tfactorforbanks * Teq for crosssec morph
	Tfac = Bi{teltopo}./Hi{teltopo}; %3. Tfac = 1/fraction of width that banks (~H) take
   Tw{teltopo} = Tfac.* Hi{teltopo}.*Bi{teltopo}./ (qi{teltopo}.*annual);
   Bequi = aB .* Q(teltopo)^bB .* D(teltopo).^bD .* M(teltopo).^bM;
   Bi{teltopo} = Bi{teltopo} + (Bequi-Bi{teltopo})./Tw{teltopo};
%   B = repmat(Bequi,size(Bprevious)); %immediate adaptation
   %Exner, conserve width mass
   etai{teltopo} = [etaitemp{teltopo}] - (1/(1-lamp)).*(dQdx{teltopo}./Bi{teltopo}).*If*annual*dt...
      + (Bi{teltopo}-Bprevious{teltopo}).*Hi{teltopo}./Bi{teltopo};
   etai{teltopo}(itect{teltopo}) = etai{teltopo}(itect{teltopo}) + etatect(teltopo)*dt;
   %hold on; plot(telm,etai(1)-etaitemp(1),'+')
   Bprevious{teltopo} = Bi{teltopo};
   etaitemp{teltopo} = etai{teltopo};
end
