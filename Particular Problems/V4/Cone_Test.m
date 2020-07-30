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

CtoS = @(t,V) ulVol+(olVol-ulVol)/(2*delta*(1-t)+(olVol-ulVol)*t)*(V-(v0-delta)-(ulVol-(v0-delta))*t);
StoC = @(t,V) (V-ulVol)*(2*delta*(1-t)+(olVol-ulVol)*t)/(olVol-ulVol)+(v0-delta)+(ulVol-(v0-delta))*t;

figure(1);
xlabel('t');
ylabel('$\hat{V}$','Interpreter','latex');
title('Cone[circles] and C(S(C))[diamonds]');
hold on;
grid on;
figure(2);
xlabel('t');
ylabel('$\hat{V}$','Interpreter','latex');
title('Cone to Square');
hold on;
grid on;

for i = 1:length(t)
    
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
        scatter(t(i),CtoS(t(i),vecVols(j)),'filled');
        Cone2Square(i,j) = CtoS(t(i),vecVols(j));
    end
    
    set(0,'CurrentFigure',1);
    for j = 1:length(vecVols)
        scatter(t(i),StoC(t(i),CtoS(t(i),vecVols(j))),100,'d');
        Cone2Square2Cone(i,j) = StoC(t(i),CtoS(t(i),vecVols(j)));
    end
    
end

set(0,'CurrentFigure',1);
saveas(gcf,'Cone_Plot','epsc');
set(0,'CurrentFigure',2);
saveas(gcf,'Square_Plot','epsc');