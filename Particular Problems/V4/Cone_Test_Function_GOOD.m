close all
clear all
clc

T = 1;
dt = 0.1;
t = 0:dt:T;
olVol = 1;
ulVol = 0.7;
v0 = 0.8;
delta = 0.01;

olFV = @(t) (v0+delta)+(olVol-v0-delta)*t/T;
ulFV = @(t) (v0-delta)+(ulVol-v0+delta)*t/T;

divitions = 20;

figure(1);
xlabel('Time');
ylabel('Normalized State');
title('State-Space: Cone');
hold on;
grid minor;
ylim([0.7 1]);
box on;
figure(2);
xlabel('Time');
ylabel('Normalized State');
title('State-Space: Square');
ylim([0.7 1]);
hold on;
grid minor;
box on;

for i = 1:length(t)
    
    if i < 4
        divitions = 4;
    elseif i < 9
        divitions = 8;
    else
        divitions = 16;
    end
    
    len = olFV(t(i)) - ulFV(t(i));
    divs = len/divitions;
    vecVols = ulFV(t(i)):divs:olFV(t(i));
    
    set(0,'CurrentFigure',1);
    for j = 1:length(vecVols)
        scatter(t(i),vecVols(j),'filled');
        Cone(i,j) = vecVols(j);
    end
    
    set(0,'CurrentFigure',2);
    for j = 1:length(vecVols)
        [a,b] = CtoS(ulVol,olVol,delta,v0,t(i),vecVols(j));
        scatter(a,b,'filled');
        Cone2Square(i,j) = b;
    end
    
    set(0,'CurrentFigure',1);
    for j = 1:length(vecVols)
        [a,b] = CtoS(ulVol,olVol,delta,v0,t(i),vecVols(j));
        [c,d] = StoC(ulVol,olVol,delta,v0,a,b);
        %scatter(c,d,100,'c');
        Cone2Square2Cone(i,j) = d;
    end
    
end

set(0,'CurrentFigure',1);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
saveas(gcf,'Cone_Plot','epsc');
set(0,'CurrentFigure',2);
set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(gcf,'Square_Plot','epsc');