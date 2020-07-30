close all
clear all
clc

%% Bonete:

Level = [16:0.25:28];
Flow = [654 644 634 624 692 662 652 642 703 692 681 670 699 688 679 669 696 687 ...
    677 668 706 697 688 679 692 684 676 668 691 683 676 666 677 670 664 656 667 661 655 649 660 654 648 642 648 642 636 620 624];

figure
plot(Level,Flow)
xlabel('Level (m)')
ylabel('Turbine Flow (m^3/s)')
title('Maximum turbine flow Vs Level of Bonete')
xlim([min(Level) max(Level)])
grid on

p = polyfit(Level,Flow,2);
hold on
y = polyval(p,Level);
plot(Level,y)
legend('Real data','Second order approximation')

E(1,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(1,2) = max(abs(Flow-y));
E(1,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(1,4) = max(abs(Flow-y)./Flow);

P{1} = p;

saveas(gcf,'Bon_2','epsc')

%% Baygorria:

Level = [8.5:0.5:17.5];
Flow = [744 732 713 711 688 836 831 825 834 857 870 901 866 820 786 723 687 656 630];

figure
plot(Level,Flow);
xlabel('Level (m)')
ylabel('Turbine Flow (m^3/s)')
title('Maximum turbine flow Vs Level of Baygorria')
xlim([min(Level) max(Level)])
grid on

p = polyfit(Level,Flow,2);
hold on
y = polyval(p,Level);
plot(Level,y)
legend('Real data','Second order approximation')

E(2,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(2,2) = max(abs(Flow-y));
E(2,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(2,4) = max(abs(Flow-y)./Flow);

P{2} = p;

saveas(gcf,'Bay_2','epsc')

%% Palmar:

Level = [16:31];
Flow = [1362 1422 1519 1520 1522 1640 1707 1703 1632 1542 1466 1400 1344 1293 1247 1206];

figure
plot(Level,Flow);
xlabel('Level (m)')
ylabel('Turbine Flow (m^3/s)')
title('Maximum turbine flow Vs Level of Palmar')
xlim([min(Level) max(Level)])
grid on

p = polyfit(Level,Flow,2);
hold on
y = polyval(p,Level);
plot(Level,y)
legend('Real data','Second order approximation')

E(3,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(3,2) = max(abs(Flow-y));
E(3,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(3,4) = max(abs(Flow-y)./Flow);

P{3} = p;

saveas(gcf,'Pal_2','epsc')

%% Bonete:

Level = [16:0.25:28];
Flow = [654 644 634 624 692 662 652 642 703 692 681 670 699 688 679 669 696 687 ...
    677 668 706 697 688 679 692 684 676 668 691 683 676 666 677 670 664 656 667 661 655 649 660 654 648 642 648 642 636 620 624];

figure
plot(Level,Flow)
xlabel('Level (m)')
ylabel('Turbine Flow (m^3/s)')
title('Maximum turbine flow Vs Level of Bonete')
xlim([min(Level) max(Level)])
grid on

p = polyfit(Level(1:21),Flow(1:21),1);
P{4} = p;
hold on
y1 = polyval(p,Level(1:21));
p = polyfit(Level(22:end),Flow(22:end),1);
P{5} = p;
hold on
y2 = polyval(p,Level(22:end));
y = [y1,y2];
plot(Level,y)
legend('Real data','Segmented first order approximation')

E(4,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(4,2) = max(abs(Flow-y));
E(4,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(4,4) = max(abs(Flow-y)./Flow);

saveas(gcf,'Bon_1','epsc')

%% Baygorria:

Level = [8.5:0.5:17.5];
Flow = [744 732 713 711 688 836 831 825 834 857 870 901 866 820 786 723 687 656 630];

figure
plot(Level,Flow);
xlabel('Level (m)')
ylabel('Turbine Flow (m^3/s)')
title('Maximum turbine flow Vs Level of Baygorria')
xlim([min(Level) max(Level)])
grid on

p = polyfit(Level(1:12),Flow(1:12),1);
P{6} = p;
hold on
y1 = polyval(p,Level(1:12));
p = polyfit(Level(13:end),Flow(13:end),1);
P{7} = p;
hold on
y2 = polyval(p,Level(13:end));
y = [y1,y2];
plot(Level,y)
legend('Real data','Segmented first order approximation')

E(5,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(5,2) = max(abs(Flow-y));
E(5,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(5,4) = max(abs(Flow-y)./Flow);

saveas(gcf,'Bay_1','epsc')

%% Palmar:

Level = [16:31];
Flow = [1362 1422 1519 1520 1522 1640 1707 1703 1632 1542 1466 1400 1344 1293 1247 1206];

figure
plot(Level,Flow);
xlabel('Level (m)')
ylabel('Turbine Flow (m^3/s)')
title('Maximum turbine flow Vs Level of Palmar')
xlim([min(Level) max(Level)])
grid on

p = polyfit(Level(1:8),Flow(1:8),1);
P{8} = p;
hold on
y1 = polyval(p,Level(1:8));
p = polyfit(Level(9:end),Flow(9:end),1);
P{9} = p;
hold on
y2 = polyval(p,Level(9:end));
y = [y1,y2];
plot(Level,y)
legend('Real data','Segmented first order approximation')

E(6,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(6,2) = max(abs(Flow-y));
E(6,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(6,4) = max(abs(Flow-y)./Flow);

saveas(gcf,'Pal_1','epsc')

%% Bonete:

Level = [16:0.25:28];
Flow = [654 644 634 624 692 662 652 642 703 692 681 670 699 688 679 669 696 687 ...
    677 668 706 697 688 679 692 684 676 668 691 683 676 666 677 670 664 656 667 661 655 649 660 654 648 642 648 642 636 620 624] ;

figure
plot(Level,Flow)
xlabel('Level (m)')
ylabel('Turbine Flow (m$^3$/s)','Interpreter','latex')
title('Maximum turbine flow Vs Level of Bonete')
xlim([min(Level) max(Level)])
grid on

p = polyfix(Level(1:20),Flow(1:20),Level(21),max(Flow),1);
P{10} = p;
hold on
y1 = polyval(p,Level(1:21));
p = polyfix(Level(22:end),Flow(22:end),Level(21),max(Flow),1);
P{11} = p;
hold on
y2 = polyval(p,Level(22:end));
y = [y1,y2];
plot(Level,y)
legend('Real data','Piecewise linear approximation','location','southeast')

E(4,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(4,2) = max(abs(Flow-y));
E(4,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(4,4) = max(abs(Flow-y)./Flow);

saveas(gcf,'Bon_15','epsc')

%% Baygorria:

Level = [8.5:0.5:17.5];
Flow = [744 732 713 711 688 836 831 825 834 857 870 901 866 820 786 723 687 656 630];
% Flow = Flow - max(Flow) + 828; % We do this modification to no overload the dam.

figure
plot(Level,Flow);
xlabel('Level (m)')
ylabel('Turbine Flow (m$^3$/s)','Interpreter','latex')
title('Maximum turbine flow Vs Level of Baygorria')
xlim([min(Level) max(Level)])
grid on

p = polyfix(Level(1:11),Flow(1:11),Level(12),max(Flow),1);
P{12} = p;
hold on
y1 = polyval(p,Level(1:12));
p = polyfix(Level(13:end),Flow(13:end),Level(12),max(Flow),1);
P{13} = p;
hold on
y2 = polyval(p,Level(13:end));
y = [y1,y2];
plot(Level,y)
legend('Real data','Piecewise linear approximation','location','southeast')

E(5,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(5,2) = max(abs(Flow-y));
E(5,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(5,4) = max(abs(Flow-y)./Flow);

saveas(gcf,'Bay_15','epsc')

%% Palmar:

Level = [16:31];
Flow = [1362 1422 1519 1520 1522 1640 1703 1707 1632 1542 1466 1400 1344 1293 1247 1206];
% Flow = Flow - max(Flow) + 1373; % We do this modification to no overload the dam.

figure
plot(Level,Flow);
xlabel('Level (m)')
ylabel('Turbine Flow (m$^3$/s)','Interpreter','latex')
title('Maximum turbine flow Vs Level of Palmar')
xlim([min(Level) max(Level)])
grid on

p = polyfix(Level(1:7),Flow(1:7),Level(8),max(Flow),1);
P{14} = p;
hold on
y1 = polyval(p,Level(1:8));
p = polyfix(Level(9:end),Flow(9:end),Level(8),max(Flow),1);
P{15} = p;
hold on
y2 = polyval(p,Level(9:end));
y = [y1,y2];
plot(Level,y)
legend('Real data','Piecewise linear approximation','location','southeast')

E(6,1) = sqrt(sum(abs(Flow-y).^2)/length(Flow));
E(6,2) = max(abs(Flow-y));
E(6,3) = sqrt(sum(abs((Flow-y)./Flow).^2)/length(Flow));
E(6,4) = max(abs(Flow-y)./Flow);

saveas(gcf,'Pal_15','epsc')