function [Uini] = Cost_Ini(WeUseLambda,WithWind,SaveFig,SaveU,ShowFig,ComputePlots,Use_Fmincon,expansion,Lag,Hat_Lag)

global Demand VG1 VG4

NT = 4*2^expansion; % Discretizations in time.
NV = 4*2^expansion; % Discretizations in water.
NW = 8; % Discretizations in wind.

% ==================== ShowFig ====================>>>
if ShowFig == 1
    set(groot,'defaultFigureVisible','on');
elseif ShowFig == 0
    set(groot,'defaultFigureVisible','off');
else
    disp('Choose ShowFig between 0 or 1.');
    return;
end

% ==================== Screensize ====================>>>
Screensize = get(0,'screensize');
vert = Screensize(4)/4;
horiz = Screensize(3)/5;
fig = {};

% ==================== ComputePlots ====================>>>
if ComputePlots == 1
    for k=1:5
    set(groot,'defaultfigureposition',[horiz*(k-1) 100 horiz vert])
    fig{k} = figure(k);
    end
    for k=6:10
    set(groot,'defaultfigureposition',[horiz*(k-6) 100+vert horiz vert])
    fig{k} = figure(k);
    end
    for k=11:15
    set(groot,'defaultfigureposition',[horiz*(k-11) 100+2*vert horiz vert])
    fig{k} = figure(k);
    end
    for k=16:20
    set(groot,'defaultfigureposition',[horiz*(k-16) 100 horiz vert])
    fig{k} = figure(k);
    end
    for k=21:21
    set(groot,'defaultfigureposition',[horiz*(k-21) 100+vert horiz vert])
    fig{k} = figure(k);
    end
end

% ==================== Time ====================>>>
absTime = 0:1/24:1; % Is needed for interpolation.
Tmax = 1; % Simulation final time.
Tmin = 0; % Simulation initial time.
dt = (Tmax-Tmin)/NT;
time = Tmin:dt:Tmax; % Simulation time.
realIniTime = 0; % 0 hours.
realFinTime = 24; % 24 hours.
RealTime = realIniTime:(realFinTime-realIniTime)/(length(time)-1):realFinTime;

% ==================== Parameters (Demand, Water and Fuel) ====================>>
FlowMax4 = 4.2*10^3; % In m^3/s.
FlowMax1 = 640; % In m^3/s.
FuelMax = 2*10^6; % In kW.
VolMax4 = 5*10^9; % In m^3.
VolMax1 = 9.5*10^9; % In m^3.
TMax = 24*3600; % In s.
DMax = 2*10^6; % In kW.
Vini4 = 0.9; % Initial condition between 0.6 and 1.
Vini1 = 0.8; % Initial condition between 0.2 and 1.
Z41n = 1/3; % Delay between Bonete and Salto Grande.
Z41 = floor(length(time)*Z41n); % Discrete delay.
sigH1 = 0;
sigH4 = 0;
load('D01012017.mat'); % Demand from 01/01/2017.
IT4 = @(t) FlowMax4/2; % In m^3/s.
IT1 = @(t) FlowMax1/2; % In m^3/s.
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
D = DMax*(0.5-(0.5)*sin(2.5*pi*time)); % In kW. (THIS (1))
% D = DMax*(0.5)*time./time; % In kW. (THIS (2))
% D = D01012017'*1000*0.8; % In kW. (THIS (3))
% D = interp1(absTime,D,time); % Interpolation. (THIS (3))

% ==================== Coding ====================>>
E = 1; % To evaluate possible errors in the Lagrange Multiplayers.
E2 = 1; % To evaluate possible errors in the optimization.
Scenario = 0; % To use during our optimization.
NonLinCounter = 0; % To use with fmincon.
NoSolution = 1; % To use with fmincon.

% ==================== Parameters (Wind) ====================>>
load('Wind_Data.mat');
GM = 1383; % MW. Maximum wind power.
Gk = GM*10^3; % kW. Maximum wind power.
ExW = Wind_Data(1,:)'; % Wind in expected value.
d_ExW = Wind_Data(2,:)'; % Derivative of [Wind in expect value].
SigW = Wind_Data(3,:)'; % Standard deviation of the wind.
d_SigW = Wind_Data(4,:)'; % Derivative of [Standard deviation of the wind].
ExW = interp1(absTime,ExW,time); % Interpolation.
d_ExW = interp1(absTime,d_ExW,time); % Interpolation.
d_pW = d_ExW; % Interpolated derivative of the forecast (for our ODE, forecast = expected value).
SigW = interp1(absTime,SigW,time); % Interpolation.
d_SigW = interp1(absTime,d_SigW,time); % Interpolation.
theta_W0 = 0.2; % W/s -> 2.0000
theta_W1 = 0.03;
alpha_W = 0.14;
theta_W = @(t) theta_W0*exp(-theta_W1*t);
pW = ExW; % Interpolated forecast (for our ODE, forecast = expected value).
MinW = 0.3; % Maximum wind (it is between 0 and 1).
MaxW = 0.7; % Minimum wind (it is between 0 and 1).
dw = (MaxW-MinW)/NW;
Y_hat = MinW:dw:MaxW; % Discretized wind (it is between 0 and 1).
for i=1:length(Y_hat)
    for j=1:length(D)
        Y(j,i) = Y_hat(i); % Y(:,i) is the i-th line of wind.
        DW(j,i) = max(0,D(j)-Y(j,i)*Gk); % DW(:,i) is the i-th line of demand.
    end
end
if WithWind == 1
    WindLow = 1;
    WindHigh = length(Y_hat);
elseif WithWind == 0
    DW(:,floor((NW)/2)+1) = D(:);
    alpha_W = 0;
    theta_W = @(t) 0;
    WindLow = floor((NW)/2)+1;
    WindHigh = floor((NW)/2)+1;
else
    disp('Choose WithWind between 0 or 1.');
    return;
end

% ==================== Plot Demand (plot 6) ====================>>
if ComputePlots == 1
    set(0,'CurrentFigure',6);
    hold on
    for i=1:length(Y_hat)
        if WithWind == 1
            plot(time,DW(:,i));
        elseif (WithWind == 0)
            if i == floor((NW)/2)+1
                plot(time,DW(:,i));
            end
        elseif WithWind ~= 0
            disp('Choose WithWind between 0 or 1.');
            return;
        end
    end
    title('Effective Demand');
    xlabel('Time');
    ylabel('kW');
    grid on;
    pause(0.001);
elseif ComputePlots ~= 0
    disp('Choose ComputePlots between 0 or 1.');
    return;
end

% ==================== Matrices ====================>>
u = zeros(length(time),length(V1),length(V4),length(Y_hat)); % 3 wind lines.
du1val = zeros(length(time),length(V1),length(V4),length(Y_hat));
du4val = zeros(length(time),length(V1),length(V4),length(Y_hat));
duWval = zeros(length(time),length(V1),length(V4),length(Y_hat));
Ph4 = zeros(length(time),length(V1),length(V4));
Ph1 = zeros(length(time),length(V1),length(V4));
Pf = zeros(length(time),length(V1),length(V4));
PhiT4 = zeros(length(time),length(V1),length(V4),length(Y_hat));
VPhiT4 = zeros(length(time),length(V1),length(V4),length(Y_hat));
PhiT1 = zeros(length(time),length(V1),length(V4),length(Y_hat));
PhiF = zeros(length(time),length(V1),length(V4),length(Y_hat));
Ph4_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
Ph1_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
Pf_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
HAM_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
HAM_Plot_HF = zeros(length(time),length(V1),length(V4),length(Y_hat));
HAM_Plot_W = zeros(length(time),length(V1),length(V4),length(Y_hat));
uw = zeros(length(V1),length(V4),length(Y_hat));
uww = zeros(length(V1),length(V4),length(Y_hat));

% ==================== More parameters (Water and Fuel) ====================>>
H4 = @(V) (-19.8)*(V).^2 + (51.5)*(V) + 3.79; % Height of the water as a function of the normal volume (dam 4).
H1 = @(V) (-3.74)*(V).^2 + (16.7)*(V) + 67.7; % Height of the water as a function of the normal volume (dam 1).
cd4 = 5.3;
d4 = @(tur) cd4*tur; % Variation in the high due of the turbined flow (dam 4).
cd1 = 1.0;
d1 = @(tur) cd1*tur; % Variation in the high due of the turbined flow (dam 1).
eta = 8.72; % Efficiency.
h04 = 5.1; % Level after Salto Grande.
h01 = 54; % Level after Bonete (this will be later the level of Baygorria).
Khe4 = 0.3e-6; % In M USD/MWh.
Khe1 = 0.3e-6; % In M USD/MWh.
Kf = 130e-6; % In M USD/MWh.
Kh4 = Khe4*eta*(H4(Vini4)-d4(1)-h04)/(3.6*10^6); % Water's value of Salto Grande.
Kh1 = Khe1*eta*(H1(Vini1)-d1(1)-h01)/(3.6*10^6); % Water's value of Bonete.
hatKf = Kf*FuelMax/(3.6*10^6);
hatKh4 = @(uv) FlowMax4*(Kh4-uv/VolMax4);
hatKh1 = @(uv) FlowMax1*(Kh1-uv/VolMax1);

% ==================== Display costs ====================>>
disp(['Bonete in Cost function: ',num2str(Kh1*TMax*FlowMax1)]);
disp(['Salto Grande in Cost function: ',num2str(Kh4*TMax*FlowMax4)]);
disp(['Fuel Station in Cost function: ',num2str(Kf*TMax*FuelMax/(3.6*10^6))]);
disp(['Bonete power cost: ',num2str(Khe1),' MUSD/MWh']);
disp(['Salto Grande power cost: ',num2str(Khe4),' MUSD/MWh']);
disp(['Fuel Station power cost: ',num2str(1000*Kf*TMax*FuelMax/(3.6*10^6)/FuelMax),' MUSD/MWh']);

% ==================== Hamiltonian ====================>>
HamToMin = @(uv4,uv1,uvv4,uvv1,uw,uww,phit4,phit1,phif,Vphit4,t,Y_hat) (uv4*TMax/VolMax4*FlowMax1+Lag(t))*Vphit4*WeUseLambda ... % Virtual control.
    + TMax*hatKh4(uv4)*phit4 ... % Phi_4 component.
    + (TMax*hatKh1(uv1)-(Hat_Lag(t+Z41))*WeUseLambda)*phit1 ... % Phi_1 component.
    + TMax*hatKf*phif ... % Fuel component.
    + TMax*IT4(t)*uv4/VolMax4 + TMax*IT1(t)*uv1/VolMax1 ... % Inflow of Water.
    + (sigH1^2)/2*uvv1 + (sigH4^2)/2*uvv4 ... % Diffusion in Water.
    + (d_pW(t)+theta_W(RealTime(t))*(pW(t)-Y_hat))*uw + theta_W(RealTime(t))*alpha_W*Y_hat*(1-Y_hat)*uww; % Wind.

HamToMin_HF = @(uv4,uv1,uvv4,uvv1,phit4,phit1,phif,Vphit4,t) (uv4*TMax/VolMax4*FlowMax1+Lag(t))*Vphit4*WeUseLambda ... % Virtual control.
    + TMax*hatKh4(uv4)*phit4 ... % Phi_4 component.
    + (TMax*hatKh1(uv1)-(Hat_Lag(t+Z41))*WeUseLambda)*phit1 ... % Phi_1 component.
    + TMax*hatKf*phif ... % Fuel component.
    + TMax*IT4(t)*uv4/VolMax4 + TMax*IT1(t)*uv1/VolMax1 ... % Inflow of Water.
    + (sigH1^2)/2*uvv1 + (sigH4^2)/2*uvv4; % Diffusion in Water.
    
HamToMin_W = @(uw,uww,t,Y_hat) (d_pW(t)+theta_W(RealTime(t))*(pW(t)-Y_hat))*uw + theta_W(RealTime(t))*alpha_W*Y_hat*(1-Y_hat)*uww; % Wind.

% ==================== Control <==> Power ====================>>
K1_4 = @(V) eta*FlowMax4*(H4(Vini4)-h04); % READ: Minimum cost condition.
K1_1 = @(V) eta*FlowMax1*(H1(V)-h01); % READ: Minimum cost condition.
K2_4 = -eta*cd4*FlowMax4; % READ: Minimum cost condition.
K2_1 = -eta*cd1*FlowMax1; % READ: Minimum cost condition.
olK2_4 = -K2_4; % READ: Minimum cost condition.
olK2_1 = -K2_1; % READ: Minimum cost condition.
phit4 = @(V,Ph) (K1_4(V)-sqrt(K1_4(V)^2-4*olK2_4*Ph))/(2*olK2_4); % Normalized control given the power and the volume (dam 4).
phit1 = @(V,Ph) (K1_1(V)-sqrt(K1_1(V)^2-4*olK2_1*Ph))/(2*olK2_1); % Normalized control given the power and the volume (dam 4).
PowH4 = @(phihatH4,V4) eta*FlowMax4*phihatH4*(H4(Vini4)-d4(phihatH4)-h04); % Power given the control and the volume (dam 4).
PowH1 = @(phihatH1,V1) eta*FlowMax1*phihatH1*(H1(V1)-d1(phihatH1)-h01); % Power given the control and the volume (dam 1).
PowF = @(phihatF) FuelMax*phihatF; % Power given the control.
olPh4 = @(V) eta*FlowMax4*(H4(Vini4)-d4(1)-h04); % Maximum power given the volume (dam 4).
olPh1 = @(V) eta*FlowMax1*(H1(V)-d1(1)-h01); % Maximum power given the volume (dam 1).

% ==================== Lagrange Multiplayers of the cost function ====================>>
A1 = @(V1)eta*FlowMax1*(H1(V1)-h01);
A4 = @(V4)eta*FlowMax4*(H4(Vini4)-h04);
B1 = eta*cd1*FlowMax1;
B4 = eta*cd4*FlowMax4;
C1 = @(uv1,phi1op,V1) hatKh1(uv1)/(A1(V1)-2*B1*phi1op);
C4 = @(uv4,phi4op,V4) hatKh4(uv4)/(A4(V4)-2*B4*phi4op);
CF = hatKf/FuelMax;

% ==================== Our optimization ====================>>
c_op = hatKf/FuelMax;
d1_op = @(uv1) hatKh1(uv1)/(2*olK2_1*sqrt(4*olK2_1));
d4_op = @(uv4) hatKh4(uv4)/(2*olK2_4*sqrt(4*olK2_4));
e1_op = @(V1) (K1_1(V1)^2)/(4*olK2_1);
e4_op = @(V4) (K1_4(V4)^2)/(4*olK2_4);
g4_op = @(uv1,uv4,V1,V4)  e4_op(V4)-e1_op(V1)*(d4_op(uv4)/d1_op(uv1))^2;
h4_op = @(uv1,uv4) (d4_op(uv4)/d1_op(uv1))^2;

% ==================== Simulation ====================>>
for t=length(time):-1:1
    
    disp(['Complete: ',num2str((length(time)-t)*100/length(time)),'%.']);
    
    uv4 = zeros(length(V1),length(V4),length(Y_hat));
    uv1 = zeros(length(V1),length(V4),length(Y_hat));
    if sigH1 == 0
        uvv4 = zeros(length(V1),length(V4),length(Y_hat));
    end
    if sigH4 == 0
        uvv1 = zeros(length(V1),length(V4),length(Y_hat));
    end
    if WithWind == 1
        uw = zeros(length(V1),length(V4),length(Y_hat));
        uww = zeros(length(V1),length(V4),length(Y_hat));
    end
    
    for k=WindLow:WindHigh
        if WithWind == 1
            disp(['Wind Case: ',num2str(k),' out of ',num2str(length(Y_hat))]);
        end
        for i=1:length(V1)
            for j=1:length(V4)

                if (i~=1) && (j~=1) && (i~=length(V1)) && (j~=length(V4)) && (k~=1) && (k~=length(Y_hat)) && (t~=length(time))
                    
                    if WithWind == 1
                        uww(i,j,k) = (u(t+1,i,j,k+1)-2*u(t+1,i,j,k)+u(t+1,i,j,k-1))/(dw^2);
                    end
                    
                    if IT1(t+1)-FlowMax1*PhiT1(t+1,i,j,k)>=0
                        uv1(i,j,k) = (u(t+1,i+1,j,k)-u(t+1,i,j,k))/(dV1);
                        du1val(t,i,j,k) = 1;
                    else
                        uv1(i,j,k) = (u(t+1,i,j,k)-u(t+1,i-1,j,k))/(dV1);
                        du1val(t,i,j,k) = -1;
                    end
                    
                    if IT4(t+1)+FlowMax1*VPhiT4(t+1,i,j,k)-FlowMax4*PhiT4(t+1,i,j,k)>=0
                        uv4(i,j,k) = (u(t+1,i,j+1,k)-u(t+1,i,j,k))/(dV4);
                        du4val(t,i,j,k) = 1;
                    else
                        uv4(i,j,k) = (u(t+1,i,j,k)-u(t+1,i,j-1,k))/(dV4);
                        du4val(t,i,j,k) = -1;
                    end
                    
                    if WithWind == 1
                        if d_pW(t+1)+theta_W(RealTime(t+1))*(pW(t+1)-Y_hat(k))>0
                            uw(i,j,k) = (u(t+1,i,j,k+1)-u(t+1,i,j,k))/(dw);
                            duWval(t,i,j,k) = 1;
                        else
                            uw(i,j,k) = (u(t+1,i,j,k)-u(t+1,i,j,k-1))/(dw);
                            duWval(t,i,j,k) = -1;
                        end
                    end
                    
                elseif t~=length(time)

                    % dV1 I
                    if i==1
                        if IT1(t+1)-FlowMax1*PhiT1(t+1,i+1,j,k)>=0
                            uv1(i,j,k) = (u(t+1,i+2,j,k)-u(t+1,i+1,j,k))/(dV1);
                            du1val(t,i,j,k) = 1;
                        else
                            uv1(i,j,k) = (u(t+1,i+1,j,k)-u(t+1,i,j,k))/(dV1);
                            du1val(t,i,j,k) = -1;
                        end
                    end
                    
                    if i==length(V1)
                        if IT1(t+1)-FlowMax1*PhiT1(t+1,i-1,j,k)>=0
                            uv1(i,j,k) = (u(t+1,i,j,k)-u(t+1,i-1,j,k))/(dV1);
                            du1val(t,i,j,k) = 1;
                        else
                            uv1(i,j,k) = (u(t+1,i-1,j,k)-u(t+1,i-2,j,k))/(dV1);
                            du1val(t,i,j,k) = -1;
                        end
                    end
                    
                    if (j==1 || j==length(V4) || k==1 || k==length(Y_hat)) && (i~=1) && (i~=length(V1))
                        if IT1(t+1)-FlowMax1*PhiT1(t+1,i,j,k)>=0
                            uv1(i,j,k) = (u(t+1,i+1,j,k)-u(t+1,i,j,k))/(dV1);
                            du1val(t,i,j,k) = 1;
                        else
                            uv1(i,j,k) = (u(t+1,i,j,k)-u(t+1,i-1,j,k))/(dV1);
                            du1val(t,i,j,k) = -1;
                        end
                    end
                    % dV1 F
                    
                    % dV4 I
                    if j==1
                        if IT4(t+1)+FlowMax1*VPhiT4(t+1,i,j+1,k)-FlowMax4*PhiT4(t+1,i,j+1,k)>=0
                            uv4(i,j,k) = (u(t+1,i,j+2,k)-u(t+1,i,j+1,k))/(dV4);
                            du4val(t,i,j,k) = 1;
                        else
                            uv4(i,j,k) = (u(t+1,i,j+1,k)-u(t+1,i,j,k))/(dV4);
                            du4val(t,i,j,k) = -1;
                        end
                    end
                    
                    if j==length(V4)
                        if IT4(t+1)+FlowMax1*VPhiT4(t+1,i,j-1,k)-FlowMax4*PhiT4(t+1,i,j-1,k)>=0
                            uv4(i,j,k) = (u(t+1,i,j,k)-u(t+1,i,j-1,k))/(dV4);
                            du4val(t,i,j,k) = 1;
                        else
                            uv4(i,j,k) = (u(t+1,i,j-1,k)-u(t+1,i,j-2,k))/(dV4);
                            du4val(t,i,j,k) = -1;
                        end
                    end
                    
                    if (i==1 || i==length(V1) || k==1 || k==length(Y_hat)) && j~=1 && j~=length(V4)
                        if IT4(t+1)+FlowMax1*VPhiT4(t+1,i,j,k)-FlowMax4*PhiT4(t+1,i,j,k)>=0
                            uv4(i,j,k) = (u(t+1,i,j+1,k)-u(t+1,i,j,k))/(dV4);
                            du4val(t,i,j,k) = 1;
                        else
                            uv4(i,j,k) = (u(t+1,i,j,k)-u(t+1,i,j-1,k))/(dV4);
                            du4val(t,i,j,k) = -1;
                        end
                    end
                    % dV4 F
                    
                    % dW I
                    if WithWind == 1
                        if k==1
                            if d_pW(t+1)+theta_W(RealTime(t+1))*(pW(t+1)-Y_hat(k+1))>0
                                uw(i,j,k) = (u(t+1,i,j,k+2)-u(t+1,i,j,k+1))/(dw);
                                duWval(t,i,j,k) = 1;
                            else
                                uw(i,j,k) = (u(t+1,i,j,k+1)-u(t+1,i,j,k))/(dw);
                                duWval(t,i,j,k) = -1;
                            end
                        end

                        if k==length(Y_hat)
                            if d_pW(t+1)+theta_W(RealTime(t+1))*(pW(t+1)-Y_hat(k-1))>0
                                uw(i,j,k) = (u(t+1,i,j,k)-u(t+1,i,j,k-1))/(dw);
                                duWval(t,i,j,k) = 1;
                            else
                                uw(i,j,k) = (u(t+1,i,j,k-1)-u(t+1,i,j,k-2))/(dw);
                                duWval(t,i,j,k) = -1;
                            end
                        end

                        if (i==1 || i==length(V1) || j==1 || j==length(V4)) && k~=1 && k~=length(Y_hat)
                            if d_pW(t+1)+theta_W(RealTime(t+1))*(pW(t+1)-Y_hat(k))>0
                                uw(i,j,k) = (u(t+1,i,j,k+1)-u(t+1,i,j,k))/(dw);
                                duWval(t,i,j,k) = 1;
                            else
                                uw(i,j,k) = (u(t+1,i,j,k)-u(t+1,i,j,k-1))/(dw);
                                duWval(t,i,j,k) = -1;
                            end
                            uww(i,j,k) = (u(t+1,i,j,k+1)-2*u(t+1,i,j,k)+u(t+1,i,j,k-1))/(dw^2);
                        end
                    end
                    % dW F
                    
                end

                % ============================ HAMILTONIAN ============================
                
                % No restrictions:
                Done = 0;                

                Ph1(t,i,j) = e1_op(V1(i))-(d1_op(uv1(i,j,k))/(2*c_op))^2;
                Ph4(t,i,j) = e4_op(V4(j))-(d4_op(uv4(i,j,k))/(2*c_op))^2;
                Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j);
                if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                    Done = 1;
                    Scenario = 1; % For debugging.
                    if uv4(i,j,k)*TMax/VolMax4*FlowMax1+Lag(t) > 0
                        VPhiT4(t,i,j,k) = 0;
                    else
                        VPhiT4(t,i,j,k) = 1;
                    end
                end

                % One restriction:
                if not(Done)
                    CostVal = [];
                    Optimums = {};
                    Optimums{1} = 0;

                    % gD_0:
                    s = (DW(t,k)-g4_op(uv1(i,j,k),uv4(i,j,k),V1(i),V4(j)))/(1+h4_op(uv1(i,j,k),uv4(i,j,k)));
                    Ph4(t,i,j) = g4_op(uv1(i,j,k),uv4(i,j,k),V1(i),V4(j)) + h4_op(uv1(i,j,k),uv4(i,j,k))*s;
                    Ph1(t,i,j) = DW(t,k)-Ph4(t,i,j);
                    Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % gD_1:
                    s = ((DW(t,k)-FuelMax)-g4_op(uv1(i,j,k),uv4(i,j,k),V1(i),V4(j)))/(1+h4_op(uv1(i,j,k),uv4(i,j,k)));
                    Ph4(t,i,j) = g4_op(uv1(i,j,k),uv4(i,j,k),V1(i),V4(j)) + h4_op(uv1(i,j,k),uv4(i,j,k))*s;
                    Ph1(t,i,j) = (DW(t,k)-FuelMax)-Ph4(t,i,j);
                    Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    %g1_0:
                    Ph1(t,i,j) = 0;
                    Ph4(t,i,j) = e4_op(V4(j)) - (d4_op(uv4(i,j,k))/(2*c_op))^2;
                    Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    %g1_1:
                    Ph1(t,i,j) = olPh1(V1(i));
                    Ph4(t,i,j) = e4_op(V4(j)) - (d4_op(uv4(i,j,k))/(2*c_op))^2;
                    Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    %g4_0:
                    Ph4(t,i,j) = 0;
                    Ph1(t,i,j) = e1_op(V1(i)) - (d1_op(uv1(i,j,k))/(2*c_op))^2;
                    Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    %g4_1:
                    Ph4(t,i,j) = olPh4(V4(j));
                    Ph1(t,i,j) = e1_op(V1(i)) - (d1_op(uv1(i,j,k))/(2*c_op))^2;
                    Pf(t,i,j) = DW(t,k)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                end

                if length(Optimums)>1
                    Done = 1;
                    if uv4(i,j,k)*TMax/VolMax4*FlowMax1+Lag(t) > 0
                        VPhiT4(t,i,j,k) = 0;
                    else
                        VPhiT4(t,i,j,k) = 1;
                    end
                    for y=2:length(Optimums)
                        CostVal(y-1) = HamToMin(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),uw(i,j,k),uww(i,j,k),phit4(V4(j),Optimums{y}(2)),phit1(V1(i),Optimums{y}(1)),Optimums{y}(3)/FuelMax,VPhiT4(t,i,j,k),t,Y_hat(k));
                    end
                    ind_op = find(CostVal == min(CostVal));
                    Ph1(t,i,j) = Optimums{ind_op+1}(1);
                    Ph4(t,i,j) = Optimums{ind_op+1}(2);
                    Pf(t,i,j) = Optimums{ind_op+1}(3);
                    Scenario = 2; % For debugging.
                end

                % Two restriction:
                if not(Done)

                    % gD_0+g1_1:
                    Ph1(t,i,j) = olPh1(V1(i));
                    Ph4(t,i,j) = DW(t,k) - Ph1(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % gD_0+g4_1:
                    Ph4(t,i,j) = olPh4(V4(j));
                    Ph1(t,i,j) = DW(t,k) - Ph4(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % g1_0+gD_1:
                    Ph1(t,i,j) = 0;
                    Ph4(t,i,j) = (DW(t,k)-FuelMax) - Ph1(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % g1_0+g4_1:
                    Ph1(t,i,j) = 0;
                    Ph4(t,i,j) = olPh4(V4(j));
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % g4_0+gD_1:
                    Ph4(t,i,j) = 0;
                    Ph1(t,i,j) = (DW(t,k)-FuelMax) - Ph4(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % g4_0+g1_1:
                    Ph4(t,i,j) = 0;
                    Ph1(t,i,j) = olPh1(V1(i));
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % gD_1+g1_1:
                    Ph1(t,i,j) = olPh1(V1(i));
                    Ph4(t,i,j) = (DW(t,k)-FuelMax) - Ph1(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % gD_1+g4_1:
                    Ph4(t,i,j) = olPh4(V4(j));
                    Ph1(t,i,j) = (DW(t,k)-FuelMax) - Ph4(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % g1_1+g4_1:
                    Ph4(t,i,j) = olPh4(V4(j));
                    Ph1(t,i,j) = olPh1(V1(i));
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % gD_0+g1_0:
                    Ph1(t,i,j) = 0;
                    Ph4(t,i,j) = DW(t,k) - Ph1(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % gD_0+g4_0:
                    Ph4(t,i,j) = 0;
                    Ph1(t,i,j) = DW(t,k) - Ph4(t,i,j);
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                    % g1_0+g4_0:
                    Ph4(t,i,j) = 0;
                    Ph1(t,i,j) = 0;
                    Pf(t,i,j) = DW(t,k) - Ph1(t,i,j) - Ph4(t,i,j);
                    if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
                        Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
                    end

                end

                if not(Done)
                    if length(Optimums)>1
                        if uv4(i,j,k)*TMax/VolMax4*FlowMax1+Lag(t) > 0
                            VPhiT4(t,i,j,k) = 0;
                        else
                            VPhiT4(t,i,j,k) = 1;
                        end

                        for y=2:length(Optimums)
                            CostVal(y-1) = HamToMin(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),uw(i,j,k),uww(i,j,k),phit4(V4(j),Optimums{y}(2)),phit1(V1(i),Optimums{y}(1)),Optimums{y}(3)/FuelMax,VPhiT4(t,i,j,k),t,Y_hat(k));
                        end
                        ind_op = find(CostVal == min(CostVal));
                            if length(ind_op)>1
                                ind_op = ind_op(1);
                                multOP = 1;
                                disp('MOP');
                            else
                                multOP = 0;
                            end
                        Ph1(t,i,j) = Optimums{ind_op+1}(1);
                        Ph4(t,i,j) = Optimums{ind_op+1}(2);
                        Pf(t,i,j) = Optimums{ind_op+1}(3);
                        Scenario = 3; % For debugging.
                    else
                        disp('Error in our Optimization, it never finishes');
                        return;
                    end
                end
                
                if t<=Z41
                    VPhiT4(t,i,j,k) = 0.5;
                end
                
                % ============================ HAMILTONIAN ============================
                PhiT4(t,i,j,k) = phit4(V4(j),Ph4(t,i,j));
                if PhiT4(t,i,j,k)+1e-6>1
                    PhiT4(t,i,j,k) = 1;
                end
                PhiT1(t,i,j,k) = phit1(V1(i),Ph1(t,i,j));
                if PhiT1(t,i,j,k)+1e-6>1
                    PhiT1(t,i,j,k) = 1;
                end
                PhiF(t,i,j,k) = Pf(t,i,j)/FuelMax;
                
                Ph1_Plot(t,i,j,k) = Ph1(t,i,j);
                Ph4_Plot(t,i,j,k) = Ph4(t,i,j);
                Pf_Plot(t,i,j,k) = Pf(t,i,j);
                                
                if t==length(time)
                    HAM(i,j) = 0;
                    HAM_Plot(t,i,j,k) = HAM(i,j);
                    HAM_Plot_HF(t,i,j,k) = HAM(i,j);
                    HAM_Plot_W(t,i,j,k) = HAM(i,j);
                else
                    HAM(i,j) = HamToMin(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),uw(i,j,k),uww(i,j,k),PhiT4(t+1,i,j,k),PhiT1(t+1,i,j,k),PhiF(t+1,i,j,k),VPhiT4(t+1,i,j,k),t+1,Y_hat(k));
                    HAM_Plot(t,i,j,k) = HAM(i,j);
                    HAM_Plot_HF(t,i,j,k) = HamToMin_HF(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),PhiT4(t+1,i,j,k),PhiT1(t+1,i,j,k),PhiF(t+1,i,j,k),VPhiT4(t+1,i,j,k),t+1);
                    HAM_Plot_W(t,i,j,k) = HamToMin_W(uw(i,j,k),uww(i,j,k),t+1,Y_hat(k));
                end
                
                % ============================ FMINCON ============================
                if Use_Fmincon == 1
                    if t~=1% && multOP == 0
                        if NonLinCounter == 100

                            while NoSolution ~= 0
                                Demand = DW(t,k);
                                VG1 = V1(i);
                                VG4 = V4(j);

                                fun = @(x) HamToMin(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),uw(i,j,k),uww(i,j,k),x(2),x(1),x(3),VPhiT4(t,i,j,k),t,Y_hat(k));

                                A = [];
                                b = [];
                                Aeq = [];
                                beq = [];
                                lb = zeros(1,3);
                                ub = ones(1,3);
                                x0 = rand(1,3);
                                if mod(NoSolution,3) == 0
                                    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',10000,'FunctionTolerance',1.0000e-10,'MaxIterations',10000,'Algorithm','interior-point','MaxFunEvals',10^6);
                                elseif mod(NoSolution,3) == 1
                                    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',10000,'FunctionTolerance',1.0000e-10,'MaxIterations',10000,'Algorithm','active-set','MaxFunEvals',10^6);
                                else
                                    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',10000,'FunctionTolerance',1.0000e-10,'MaxIterations',10000,'Algorithm','sqp','MaxFunEvals',10^6);
                                end
                                [x,~,flagFmincon,outputFmincon] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

                                if abs(x(1)-PhiT1(t,i,j,k))>=0.01 || abs(x(2)-PhiT4(t,i,j,k))>=0.01 || abs(x(3)-PhiF(t,i,j,k))>=0.001
                                    NoSolution = NoSolution + 1;
                                    disp('The optimizations do not match!');
                                    if NoSolution == 10000
                                        disp('Maximum number of iterations reached!');
                                        disp(['Our solution (sum(P)-DW(t,w))=',num2str(Ph4(t,i,j)+Ph1(t,i,j)+Pf(t,i,j)-DW(t,k))]);
                                        disp(['Our solution H = ',num2str(HamToMin(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),uw(i,j,k),uww(i,j,k),PhiT4(t,i,j,k),PhiT1(t,i,j,k),PhiF(t,i,j,k),t,Y_hat(k),WeUseLambda,Lag,Hat_Lag))]);
                                        disp(['MATLAB solution (sum(P)-DW(t,w))=',num2str(PowH1(x(1),V1(i))+PowH4(x(2),V4(j))+PowF(x(3))-DW(t,k))]);
                                        disp(['MATLAB solution H = ',num2str(HamToMin(uv4(i,j,k),uv1(i,j,k),uvv4(i,j,k),uvv1(i,j,k),uw(i,j,k),uww(i,j,k),x(2),x(1),x(3),t,Y_hat(k),WeUseLambda,Lag,Hat_Lag))]);
                                        error2V1(E2) = abs(x(1)-PhiT1(t,i,j,k));
                                        error2V4(E2) = abs(x(2)-PhiT4(t,i,j,k));
                                        error2F(E2) = abs(x(3)-PhiF(t,i,j,k));
                                        return;
                                    end
                                else
                                    disp('The optimizations match!');
                                    error2V1(E2) = abs(x(1)-PhiT1(t,i,j,k));
                                    error2V4(E2) = abs(x(2)-PhiT4(t,i,j,k));
                                    error2F(E2) = abs(x(3)-PhiF(t,i,j,k));
                                    E2 = E2+1;
                                    NoSolution = 0;
                                end
                            end
                            NonLinCounter = 0;
                            NoSolution = 1;
                        else
                            NonLinCounter = NonLinCounter + 1;
                        end
                    end
                elseif Use_Fmincon ~= 0
                    disp('Choose Use_fmincon between 0 or 1');
                    return;
                end
                % ============================ FMINCON ============================

                PH1logic = PhiT1(t,i,j) >= 0.99999; % 1=>1.
                PH4logic = PhiT4(t,i,j) >= 0.99999; % 1=>1.
                PHFlogic = 1 - (PhiF(t,i,j) == 0); % 0=>0.

                % === LAMBDA ===>>>
                    if PH1logic==1
                        if PH4logic==1
                            if PHFlogic==0
                                Lambda(t,i,j,k) = CF;
                            elseif PHFlogic~=0
                                Lambda(t,i,j,k) = CF;
                            end
                        elseif PH4logic~=1
                            if PHFlogic==0
                                Lambda(t,i,j,k) = C4(uv4(i,j,k),PhiT4(t,i,j),V4(j));
                            elseif PHFlogic~=0
                                Lambda(t,i,j,k) = C4(uv4(i,j,k),PhiT4(t,i,j),V4(j));
                                error(E) = abs(C4(uv4(i,j,k),PhiT4(t,i,j),V4(j)) - CF)/CF;
                                E = E+1;
                            end
                        end
                    elseif PH1logic~=1
                        if PH4logic==1
                            if PHFlogic==0
                                Lambda(t,i,j,k) = C1(uv1(i,j,k),PhiT1(t,i,j),V1(i));
                            elseif PHFlogic~=0
                                Lambda(t,i,j,k) = C1(uv1(i,j,k),PhiT1(t,i,j),V1(i));
                                error(E) = abs(C1(uv1(i,j,k),PhiT1(t,i,j),V1(i)) - CF)/CF;
                                E = E+1;
                            end
                        elseif PH4logic~=1
                            if PHFlogic==0
                                Lambda(t,i,j,k) = C1(uv1(i,j,k),PhiT1(t,i,j),V1(i));
                                error(E) = C1(uv1(i,j,k),PhiT1(t,i,j),V1(i)) - C4(uv4(i,j,k),PhiT4(t,i,j),V4(j));
                                E = E+1;
                            elseif PHFlogic~=0
                                Lambda(t,i,j,k) = CF;
                                error(E) = abs((2/3)*C1(uv1(i,j,k),PhiT1(t,i,j),V1(i)) - C4(uv4(i,j,k),PhiT4(t,i,j),V4(j))/3 - CF/3)/CF;
                            end
                        end
                    end
                % === LAMBDA ===>>>

            end
        end
        if t~=length(time)
            u(t,:,:,k) = u(t+1,:,:,k) + dt*reshape(HAM,[1,length(V1),length(V4)]);
        end
        
    end

    % === PLOT ===>>>
    if ComputePlots == 1
        if WithWind == 0
            hmin = 5;
            hm = 5;
            hmax = 5;
        else
            hmin = 1;
            hm = 1;
            hmax = length(Y_hat);
        end

        set(0,'CurrentFigure',6);
        vline(time(t),'r',' ');
        title(['Demand (time=',num2str(t),')']);

        set(0,'CurrentFigure',1);
        [v1,v4] = meshgrid(V1,V4);
        clf(1);
        hold on;
        for k=hmin:hm:hmax
            s1 = surf(v1,v4,squeeze(u(t,:,:,k))');
            s1.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        zlabel('M USD');
        title(['Optimal cost function (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        box;
        view(140,20);

        set(0,'CurrentFigure',2);
        clf(2);
        hold on;
        for k=hmin:hm:hmax
            s2 = surf(v1,v4,squeeze(VPhiT4(t,:,:,k))');
            s2.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Virtual flow 4  (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',7);
        clf(7);
        hold on;
        for k=hmin:hm:hmax
            s7 = surf(v1,v4,squeeze(Lambda(t,:,:,k))');
            s7.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Lagrangian Value (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',5);
        clf(5);
        hold on;
        for k=hmin:hm:hmax
            s5 = surf(v1,v4,squeeze(PhiT1(t,:,:,k))');
            s5.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Turbined flow 1 (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',4);
        clf(4);
        hold on;
        for k=hmin:hm:hmax
            s4 = surf(v1,v4,squeeze(PhiT4(t,:,:,k))');
            s4.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Turbined flow 4 (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',3);
        clf(3);
        hold on;
        for k=hmin:hm:hmax
            s3 = surf(v1,v4,squeeze(PhiF(t,:,:,k))');
            s3.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Fuel control (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',8);
        clf(8);
        hold on;
        for k=hmin:hm:hmax
            s8 = surf(v1,v4,squeeze(duWval(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['duWval (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',10);
        clf(10);
        hold on;
        for k=hmin:hm:hmax
            s10 = surf(v1,v4,squeeze(Ph1_Plot(t,:,:,k))');
            s10.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Hydraulic Power 1 (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',9);
        clf(9);
        hold on;
        for k=hmin:hm:hmax
            s9 = surf(v1,v4,squeeze(Ph4_Plot(t,:,:,k))');
            s9.EdgeColor = 'none';
        end
        xlabel('V1');
        ylabel('V4');
        title(['Hydraulic Power 4 (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',11);
        clf(11);
        hold on;
        for k=hmin:hm:hmax
            s11 = surf(v1,v4,squeeze(uv1(:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['uv1 (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',12);
        clf(12);
        hold on;
        for k=hmin:hm:hmax
            s12 = surf(v1,v4,squeeze(uv4(:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['uv4 (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',13);
        clf(13);
        hold on;
        for k=hmin:hm:hmax
            s13 = surf(v1,v4,squeeze(HAM_Plot(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['HAMILTONIAN (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',14);
        clf(14);
        hold on;
        for k=hmin:hm:hmax
            s14 = surf(v1,v4,squeeze(HAM_Plot_HF(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['HAMILTONIAN HF (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',15);
        clf(15);
        hold on;
        for k=hmin:hm:hmax
            s15 = surf(v1,v4,squeeze(HAM_Plot_W(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['HAMILTONIAN W (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',16);
        clf(16);
        hold on;
        for k=hmin:hm:hmax
            s16 = surf(v1,v4,squeeze(du1val(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['du1val (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',17);
        clf(17);
        hold on;
        for k=hmin:hm:hmax
            s17 = surf(v1,v4,squeeze(du4val(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['du4val (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',18);
        clf(18);
        hold on;
        for k=hmin:hm:hmax
            s18 = surf(v1,v4,squeeze(duWval(t,:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['duWval (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',19);
        clf(19);
        hold on;
        for k=hmin:hm:hmax
            s19 = surf(v1,v4,squeeze(uw(:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['uw (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;

        set(0,'CurrentFigure',20);
        clf(20);
        hold on;
        for k=hmin:hm:hmax
            s20 = surf(v1,v4,squeeze(uww(:,:,k))');
        end
        xlabel('V1');
        ylabel('V4');
        title(['uww (time=',num2str(t),')']);
        xlim([Vmin1 Vmax1]);
        ylim([Vmin4 Vmax4]);
        zlim auto
        view(140,20);
        box;
        
        if t <= length(time)-Z41
            set(0,'CurrentFigure',21);
            clf(21);
            hold on;
            for k=hmin:hm:hmax
                s21 = surf(v1,v4,(squeeze(VPhiT4(t+Z41,:,:,k))'-squeeze(PhiT1(t,:,:,k))'));
            end
            xlabel('V1');
            ylabel('V4');
            title({'Difference between';['virtual and real flow (time=',num2str(t+Z41),')']});
            xlim([Vmin1 Vmax1]);
            ylim([Vmin4 Vmax4]);
            zlim([-1 1]);
            view(140,20);
            box;
            if SaveFig == 1
                saveas(gcf,['3-Difference_at_',num2str(t+Z41)],'epsc');
            end
        end

        if WithWind == 1
            if t==length(time)
                set(groot,'defaultfigureposition',[0 0 Screensize(3)/2 Screensize(4)/2]);
                fig{101} = figure(101);
                fig{102} = figure(102);
                fig{103} = figure(103);
                sV1 = floor(length(V1)/5);
                sV4 = floor(length(V4)/5);
            end

            set(0,'CurrentFigure',101);
            title(['Cost Vs Wind (time=',num2str(t),')']);
            xlabel('Normalized wind (0,1)');
            ylabel('M USD');
            grid on;
            hold on;
            plot(Y_hat,squeeze(u(t,length(V1)/2+0.5,length(V4)/2+0.5,:)),'-*');
            if SaveFig == 1
                saveas(gcf,['Cost_Wind_',num2str(t)],'epsc');
            elseif SaveFig ~= 0
                disp('Choose SaveFig between 0 or 1');
                return;
            end

            set(0,'CurrentFigure',102);
            title(['uw(w) (time=',num2str(t),')']);
            xlabel('Normalized wind (0,1)');
            ylabel('M USD');
            grid on;
            hold on;
            plot(Y_hat,squeeze(uw(length(V1)/2+0.5,length(V4)/2+0.5,:)),'-*');
            if SaveFig == 1
                saveas(gcf,['Cost_Wind_',num2str(t)],'epsc');
            elseif SaveFig ~= 0
                disp('Choose SaveFig between 0 or 1');
                return;
            end

            set(0,'CurrentFigure',103);
            title(['uww(w) (time=',num2str(t),')']);
            xlabel('Normalized wind (0,1)');
            ylabel('M USD');
            grid on;
            hold on;
            plot(Y_hat,squeeze(uww(length(V1)/2+0.5,length(V4)/2+0.5,:)),'-*');
            if SaveFig == 1
                saveas(gcf,['Cost_Wind_',num2str(t)],'epsc');
            elseif SaveFig ~= 0
                disp('Choose SaveFig between 0 or 1');
                return;
            end
        end
        pause(0.5);
    end

end

if SaveU == 1
    save(['u_V1_',num2str(length(V1)-1),'_V4_',num2str(length(V4)-1),'_W_',num2str(length(Y_hat)-1),'_T_',num2str(NT),'.mat',],'u');
    Controls = {PhiT1,PhiT4,VPhiT4,PhiF};
    save('Optimal_Controls','Controls');
elseif SaveU ~= 0
    disp('Choose SaveU between 0 or 1');
    return;
end

Uini = u(1,floor((NV)/2)+1,floor((NV)/2)+1,floor((NW)/2)+1);
return;