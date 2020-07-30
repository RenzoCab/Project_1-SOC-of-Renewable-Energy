close all;
clear all;
clc;

% ==================== DATA ====================>>>

set(groot,'defaultFigureVisible','on');
expansion = 5;
NT = 4*2^expansion; % Discretizations in time.
NV = 4*2^expansion; % Discretizations in water.
TMax = 24*3600; % In s.
VolMax4 = 5*10^9; % In m^3.
VolMax1 = 9.5*10^9; % In m^3.
Vini4 = 0.9; % Initial condition between 0.6 and 1.
Vini1 = 0.7; % Initial condition between 0.2 and 1.
FlowMax4 = 4.2*10^3; % In m^3/s.
FlowMax1 = 640; % In m^3/s.
FuelMax = 2*10^6; % In kW.
IT4 = @(t) FlowMax4/2; % In m^3/s.
IT1 = @(t) FlowMax1/2; % In m^3/s.
Tmax = 1; % Simulation final time.
Tmin = 0; % Simulation initial time.
dt = (Tmax-Tmin)/NT;
time = Tmin:dt:Tmax; % Simulation time.
Z41n = 1/3; % Delay between Bonete and Salto Grande.
Z41 = floor(length(time)*Z41n); % Discrete delay.

Khe4 = 0.3e-6; % In M USD/MWh.
Khe1 = 0.3e-6; % In M USD/MWh.
H4 = @(V) (-19.8)*(V).^2 + (51.5)*(V) + 3.79; % Height of the water as a function of the volume (dam 4).
H1 = @(V) (-3.74)*(V).^2 + (16.7)*(V) + 67.7; % Height of the water as a function of the volume (dam 1).
cd4 = 5.3;
d4 = @(tur) cd4*tur; % Variation in the high due of the turbined flow (dam 4).
cd1 = 1.0;
d1 = @(tur) cd1*tur; % Variation in the high due of the turbined flow (dam 1).
eta = 8.72; % Efficiency.
h04 = 5.1; % Level after Salto Grande.
h01 = 54; % Level after Bonete (this will be later the level of Baygorria).
DMax = 2*10^6; % In kW.
D = DMax*(0.5-(0.5)*sin(2.5*pi*time)); % In kW. (THIS (1))
PowH4 = @(phihatH4,V4) eta*FlowMax4*phihatH4*(H4(V4)-d4(phihatH4)-h04); % Power given the control and the volume (dam 4).
PowH1 = @(phihatH1,V1) eta*FlowMax1*phihatH1*(H1(V1)-d1(phihatH1)-h01); % Power given the control and the volume (dam 1).
PowF = @(phihatF) FuelMax*phihatF; % Power given the control.

Vmax4 = (Vini4*VolMax4+IT4(1)*TMax)/VolMax4;
Vmin4 = (Vini4*VolMax4-IT4(1)*TMax)/VolMax4;
Vmax1 = (Vini1*VolMax1+IT1(1)*TMax)/VolMax1;
Vmin1 = (Vini1*VolMax1-IT1(1)*TMax)/VolMax1;
DV4 = Vmax4-Vmin4;
DV1 = Vmax1-Vmin1;
dV4 = DV4/NV;
dV1 = DV1/NV;
V4 = Vmin4:dV4:Vmax4; % Discretized volume (dam 4).
V1 = Vmin1:dV1:Vmax1; % Discretized volume (dam 4).

Kh1 = Khe1*eta*(H1(Vini1)-d1(1)-h01)/(3.6*10^6); % Water's value of Bonete.
Kh4 = Khe4*eta*(H4(Vini4)-d4(1)-h04)/(3.6*10^6); % Water's value of Salto Grande.
Kf = 130e-6; % In M USD/MWh.

% ==================== NEW ====================>>>

Cf = Kf*TMax*FuelMax/(3.6e6);
Ch1 = Kh1*TMax*FlowMax1;
Ch4 = Kh4*TMax*FlowMax4;

load('Optimal_Controls');
PhiT1 = Controls{1}(:,:,:,5);
PhiT4 = Controls{2}(:,:,:,5);
VPhiT4 = Controls{3}(:,:,:,5);
PhiF = Controls{4}(:,:,:,5);

Vo1(1) = Vini1;
Vo4(1) = Vini4;
PowH1P(1) = 0;
PowH4P(1) = 0;
PowFP(1) = 0;
PhiH1(1) = 0;
VPhiH4(1) = 0;
ConH1(1) = 0;
ConH4(1) = 0;
ConF(1) = 0;

for t=1:length(time)-1
    
    [~,i1] = min(abs(Vo1(end)-V1));
    [~,i4] = min(abs(Vo4(end)-V4));
    
    PowH1P(end) = PowH1(PhiT1(t,i1,i4),Vo1(end));
    PowH4P(end) = PowH4(PhiT4(t,i1,i4),Vo4(end));
    PowFP(end) = PowF(PhiF(t,i1,i4));
    PhiH1(end) = PhiT1(t,i1,i4);
    VPhiH4(end) = VPhiT4(t,i1,i4);
    
    ConH1(end) = PhiT1(t,i1,i4);
    ConH4(end) = PhiT4(t,i1,i4);
    ConF(end) = PhiF(t,i1,i4);
    
    Vo1(end+1) = Vo1(end) + (dt*TMax/VolMax1)*(IT1(t)-PhiT1(t,i1,i4)*FlowMax1);
    Vo4(end+1) = Vo4(end) + (dt*TMax/VolMax4)*(IT4(t)+VPhiT4(t,i1,i4)*FlowMax1-PhiT4(t,i1,i4)*FlowMax4);
    
    PowH1P(end+1) = 0;
    PowH4P(end+1) = 0;
    PowFP(end+1) = 0;
    PhiH1(end+1) = 0;
    VPhiH4(end+1) = 0;
    ConH1(end+1) = 0;
    ConH4(end+1) = 0;
    ConF(end+1) = 0;
    
end

[~,i1] = min(abs(Vo1(end)-V1));
[~,i4] = min(abs(Vo4(end)-V4));
PowH1P(end) = PowH1(PhiT1(length(time),i1,i4),Vo1(end));
PowH4P(end) = PowH4(PhiT4(length(time),i1,i4),Vo4(end));
PowFP(end) = PowF(PhiF(length(time),i1,i4));
PhiH1(end) = PhiT1(length(time),i1,i4);
VPhiH4(end) = VPhiT4(length(time),i1,i4);
ConH1(end) = PhiT1(length(time),i1,i4);
ConH4(end) = PhiT4(length(time),i1,i4);
ConF(end) = PhiF(length(time),i1,i4);

t = 1:length(time);
plot(t,D/1000000,'*');
hold on;
area([PowH1P;PowH4P;PowFP]'/1000000);
legend('Demand','Power Bonete','Power SG','Power Fuel');
title('Demand and Power');
xlabel('Discretization');
grid on;
ylabel('GW');
xlim([1,length(time)]);
saveas(gcf,'Demand_and_Power','epsc');

figure;
plot(t,Vo1);
title('Volume Bonete');
xlabel('Discretization');
grid on;
ylabel('Normalized volume');
xlim([1,length(time)]);
saveas(gcf,'Volume_Bonete','epsc');

figure;
plot(t,Vo4);
title('Volume Salto Grande');
xlabel('Discretization');
ylabel('Normalized volume');
grid on;
xlim([1,length(time)]);
saveas(gcf,'Volume_Salto_Grande','epsc');

figure;
PhiH1_Sh = circshift(PhiH1,Z41);
plot(t,[zeros(1,Z41),ones(1,length(time)-Z41)].*(VPhiH4-PhiH1_Sh));
title('Virtual Flow - Real Flow');
xlabel('Discretization');
ylabel('\epsilon(t)');
xlim([1,length(time)]);
grid on;
saveas(gcf,'Virtual_Flow_Real_Flow','epsc');

figure;
plot(t,ConH1);
hold on;
plot(t,ConH4);
plot(t,ConF);
plot(t,VPhiH4,'*');
legend('Control Bonete','Control Salto Grande','Control Fuel','Virtual Control');
title('Controls');
xlabel('Discretization');
grid on;
ylabel('Normalized controls');
xlim([1,length(time)]);
saveas(gcf,'Controls','epsc');