%DO Channel belt model Network version takke44
clc
%cd 'E:\Matlabdata\ChanBmodel\outputsyst3'
%cd 'D:\Matlabdata\ChanBmodel\'
cd 'C:\Users\Maarten Kleinhans\Documents\Matlabdata\ChanBmodel\Takke44'

figure
set(gcf,'units','centimeters','position',[1 1 18 18],'papertype','A4',...
   'papertype','A4','paperunits','centimeters','paperposition',[1 1 18 18]);
magn = 0.025.*[-1 -1 1 1]; %x, y, width, height

telorder = 2;

%% A. Base case: NO gradient advantage, NO tectonics and NO bend
%cc
in44Ganga
Bifurcations(2,5) = 0;
etatect = zeros(size(etatect));
Durat = 3;
takke44
%figs44

subplot(3,3,1)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(2:3,1) = [25000 20000]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('base case','nearly symmetric'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(A) Base case')
%xlabel('time (years)')
ylabel('Relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% B. gradient advantage, NO tectonics and NO bend
%cc
in44Ganga
Sizes(3,2) = 225000;
Bifurcations(2,5) = 0;
etatect = zeros(size(etatect));
takke44
%figs44

subplot(3,3,2)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(3,2) = 240000;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'g--',...
   'linewidth',No+2-telorder)

Sizes(3,2) = 250000;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r-.',...
   'linewidth',No+2-telorder)

Sizes(3,2) = 275000;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'m:',...
   'linewidth',No+2-telorder)

h = legend('75 km adv','60','50','25'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(B) Gradient advantage')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% C. NO gradient advantage, tectonics and NO bend
%cc
in44Ganga
Bifurcations(2,5) = 0;
etatect = etatect./2;
takke44
%figs44

subplot(3,3,3)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
etatect = 2.*etatect;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'g--',...
   'linewidth',No+2-telorder)

etatect = 2.*etatect;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r-.',...
   'linewidth',No+2-telorder)


h = legend('0.025 m/yr','0.05','0.1'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(C) Subsidence')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% D. NO gradient advantage, NO tectonics and migrating bend
%cc
in44Ganga
etatect = zeros(size(etatect));
takke44
%figs44

subplot(3,3,4)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(2:3,1) = [25000 20000]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('bends','more sym'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(D) Migrating bends')
%xlabel('time (years)')
ylabel('Relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% E. gradient advantage, tectonics and NO migrating bend
%cc
in44Ganga
Sizes(3,2) = 275000;
Bifurcations(2,5) = 0;
takke44
%figs44

subplot(3,3,5)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(3,2) = 250000;
etatect = etatect./2;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'g--',...
   'linewidth',No+2-telorder)

etatect = 2.*etatect;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r-.',...
   'linewidth',No+2-telorder)

h = legend('25 km, 5 cm/yr','50, 2.5','50, 5'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(E) Grad. advan. and subs.')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% F. gradient advantage, NO tectonics and migrating bend
%cc
in44Ganga
Sizes(3,2) = 250000;
etatect = zeros(size(etatect));
takke44
%figs44

subplot(3,3,6)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(3,2) = 240000;
Bifurcations(2,5) = 0;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'g--',...
   'linewidth',No+2-telorder)

Bifurcations(2,5) = -10000;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r-.',...
   'linewidth',No+2-telorder)

h = legend('50 km','60 no bends','60'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(F) Grad. advan. and migrat. bend')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% G. NO gradient advantage, tectonics and migrating bend
%cc
in44Ganga
takke44
%figs44

subplot(3,3,7)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
etatect = 2.*etatect;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('0.05 m/yr','0.1'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(G) Subs. and migrat. bend')
xlabel('Time (years)')
ylabel('Relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% H. gradient advantage, tectonics and migrating bend
%cc
in44Ganga
Sizes(3,2) = 250000;
Bifurcations(2,5) = -1*Bifurcations(2,5);
takke44
%figs44

subplot(3,3,8)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Bifurcations(2,5) = -1*Bifurcations(2,5);
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('50 km','opposite'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(H) Preferred scenario')
xlabel('Time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% I. Width

subplot(3,3,9)
set(gca,'position',get(gca,'position')+magn,'YAxisLocation','right');

hold on
plot(rep{end}.timeser.timeax(:,1),...
  rep{end}.timeser.width1(:,(2)) ...
  ./ rep{end}.timeser.width1(:,(1)), ...
  rep{end}.timeser.timeax(:,1),...
  rep{end}.timeser.width1(:,(3)) ...
  ./ rep{end}.timeser.width1(:,(1)), ...
    'linewidth',2)
  %rep{end}.timeser.timeax(:,1),...
  %(rep{end}.timeser.width1(:,(2)) + ...
  %rep{end}.timeser.width1(:,(3)) ), ...
plot(rep{end}.timeser.timeax(:,1),...
  rep{end}.timeser.width2(:,(2)) ...
  ./ rep{end}.timeser.width2(:,(1)), ...
  rep{end}.timeser.timeax(:,1),...
  rep{end}.timeser.width2(:,(3)) ...
  ./ rep{end}.timeser.width2(:,(1)), ...
  'linewidth',1)
  %rep{end}.timeser.timeax(:,1),...
  %(rep{end}.timeser.width2(:,(2)) + ...
  %rep{end}.timeser.width2(:,(3)) ), ...
hold off
title('(I) Widths for scenario (H)')
xlabel('Time (years)')
ylabel('Relative width (-)')
axis([0 Durat 0 1])
box on


%% finish figure
pname = ['D:\wordfiles\Niladri\niladrimodelR2'];
print('-djpeg', '-r600', pname);



%% Now do Other scenarios

figure
set(gcf,'units','centimeters','position',[1 1 18 18],'papertype','A4',...
   'papertype','A4','paperunits','centimeters','paperposition',[1 1 18 18]);
magn = 0.025.*[-1 -1 1 1]; %x, y, width, height

telorder = 2;

%% A. Discharge
%cc
in44Ganga
Sizes(3,2) = 250000;
Sizes(:,1) = 60000.*[1 40/45 5/45]';
takke44
%figs44

subplot(3,3,1)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(:,1) = 30000.*[1 40/45 5/45]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

Sizes(:,1) = [45000 35000 10000]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'g-.',...
   'linewidth',No+2-telorder)

Sizes(:,1) = [45000 43000 2000]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'m-.',...
   'linewidth',No+2-telorder)

h = legend('high','low','diff asym'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(A) \Delta Discharge')
%xlabel('time (years)')
ylabel('Relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% B. Subsidence
%cc
in44Ganga
Sizes(3,2) = 250000;
etatect = etatect./2;
takke44
%figs44

subplot(3,3,2)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
etatect = 4.*etatect;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('0.025 m/yr','0.1 m/yr'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(B) \Delta Subsidence')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% C. Width
%cc
in44Ganga
Sizes(3,2) = 250000;
Bup = 2100;
takke44
%figs44

subplot(3,3,3)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Bup = 1500;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('1800 m','1500 m'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(C) \Delta Width')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% D. Gradient
%cc
in44Ganga
Sizes(3,2) = 250000;
Bifurcations(:,1) = [16 12]';
takke44
%figs44

subplot(3,3,4)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Bifurcations(:,1) = [24 18]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('4x10^{-5}','6x10^{-5}'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(D) \Delta Gradient')
%xlabel('time (years)')
ylabel('Relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% E. Flow resistance
%cc
in44Ganga
Sizes(3,2) = 250000;
kc = 0.5;
takke44
%figs44

subplot(3,3,5)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
kc = 2;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('k_s=0.5 m','k_s=2 m'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(E) \Delta Flow resistance')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% F. Particle size
%cc
in44Ganga
Sizes(3,2) = 250000;
Sizes(:,3) = [.1 .1 .1]';
takke44
%figs44

subplot(3,3,6)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
Sizes(:,3) = [.3 .3 .3]';
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('0.1 mm','0.3 mm'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(F) \Delta Particle size')
%xlabel('time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% G. Spiral flow
%cc
in44Ganga
Sizes(3,2) = 250000;
epsilon = 1;
takke44
%figs44

subplot(3,3,7)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
epsilon = 3;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('weak','strong'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(G) \Delta Bend flow')
xlabel('Time (years)')
ylabel('Relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% H. Hydraulic geometry
%cc
in44Ganga
Sizes(3,2) = 250000;
bB = 0.35;
CourantCrit = 0.5;
takke44
%figs44

subplot(3,3,8)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
bB = 0.65;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r--',...
   'linewidth',No+2-telorder)

h = legend('\beta=0.35','\beta=0.65'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(H) \Delta Hydraulic geometry')
xlabel('Time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% I. Alpha cells
%cc
in44Ganga
Sizes(3,2) = 250000;
takke44
%figs44

subplot(3,3,9)
set(gca,'position',get(gca,'position')+magn);
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),...
   'linewidth',No+2-telorder)

hold on
alw = 1;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'g--',...
   'linewidth',No+2-telorder)

alw = 3;
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'r-.',...
   'linewidth',No+2-telorder)

h = legend('pref, \alpha=2','\alpha=1','\alpha=3'); 
set(h,'box','off','fontsize',8,'Location','NorthWest');
title('(I) \Delta Length bifurcation')
xlabel('Time (years)')
%ylabel('relative discharge (-)')
axis([0 Durat 0 1]) %Sizes(1,1)


%% finish figure
pname = ['D:\wordfiles\Niladri\niladrimodelotherR2'];
print('-djpeg', '-r600', pname);



%% Test stability
figure
set(gcf,'units','centimeters','position',[1 1 9 9],'papertype','A4',...
   'papertype','A4','paperunits','centimeters','paperposition',[1 1 9 9]);

%cc
in44Ganga
Durat = 300;
switchBvar = 1;
switchperturb = 1;
Bifurcations(2,5) = 0;
etatect = zeros(size(etatect));
Sizes(2:3,1) = [Sizes(1,1)/2 Sizes(1,1)/2]';
takke44
%figs44

plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,1):Orde(telorder,2)) -Sizes(1,1)/2,...
   'linewidth',No+2-telorder)

hold on
in44Ganga
Durat = 300;
switchBvar = 1;
switchperturb = 1;
Bifurcations(2,5) = 0;
etatect = zeros(size(etatect));
Sizes(2:3,1) = [Sizes(1,1)/2 Sizes(1,1)/2]';
Zperturb = Sizes(3,3)/1000;
waterlevelprecision = 1e-6; %1e-4
dischargeprecision = 1e-4; %1e-2
takke44
plot(rep{end}.timeser.timeax(:,1),...
   rep{end}.timeser.Q(:,Orde(telorder,1):Orde(telorder,2)) -Sizes(1,1)/2,'--',...
   'linewidth',No+2-telorder)


plot([Tperturb Tperturb],[-.3 .3],'r-','linewidth',4)


% etatect = zeros(size(etatect));
% Sizes(2:3,1) = [22499.5 22500.5]';
% switchperturb = 1;
% Zperturb = 3;
% Tperturb = 50;
% Bperturb = 2;
% Nperturb = 2:250000/dx-1;
% takke44
% plot(rep{end}.timeser.timeax(:,1),...
%    rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'m:',...
%    'linewidth',No+2-telorder)
% 
% etatect = zeros(size(etatect));
% Sizes(2:3,1) = [22499.5 22500.5]';
% switchperturb = 1;
% Zperturb = -3;
% Tperturb = 50;
% Bperturb = 2;
% Nperturb = 2:250000/dx-1;
% takke44
% plot(rep{end}.timeser.timeax(:,1),...
%    rep{end}.timeser.Q(:,Orde(telorder,2)) ./Sizes(1,1),'m:',...
%    'linewidth',No+2-telorder)



h = legend('left, 1 mm','right, 1 mm','left, 0.14 mm','right, 0.14 mm',...
    'perturb time'); 
set(h,'box','off','fontsize',10,'Location','NorthWest');
title('Perturbed symmetrical system')
xlabel('Time (years)')
ylabel('Discharge difference (m^3 s^{-1})')
axis([0 Durat -inf inf]) %Sizes(1,1)

pname = ['D:\wordfiles\Niladri\niladrimodelperturbR'];
print('-djpeg', '-r600', pname);
