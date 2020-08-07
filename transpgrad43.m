%transport gradients for use in takke43

dQdx{teltopo}(1) = au*(Qi{teltopo}(1)-Qsnode(teltopo))/dx + ...
   (1-au)*(Qi{teltopo}(2)-Qi{teltopo}(1))/dx; %AMBIENT ghost node transport
for telx=2:Nx(teltopo) %transport gradient, positive means erosion
   dQdx{teltopo}(telx) = au*(Qi{teltopo}(telx)-Qi{teltopo}(telx-1))/dx + ...
      (1-au)*(Qi{teltopo}(telx+1)-Qi{teltopo}(telx))/dx;
end %transport gradient
dQdx{teltopo}(Nx(teltopo)+1) = (Qi{teltopo}(Nx(teltopo)+1)-Qi{teltopo}(Nx(teltopo)))/dx;
