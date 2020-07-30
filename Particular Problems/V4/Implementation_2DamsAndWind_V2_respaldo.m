close all
clear all
clc

global Demand VG1 VG4

set(groot,'defaultFigureVisible','off') % Show plots OFF.

E = 1; % To evaluate possible errors.
E2 = 1; % To evaluate possible errors in the optimization.
load('D01012017.mat');
load('WP01012017.mat');
load('Wind_Data.mat');
GM = 1383; % MW.
Wind_For = WP01012017/GM; % Real value of the wind.
SigW = Wind_Data(1,:)';
ExW = zeros(1,25);
d_SigW = Wind_Data(2,:)';
d_ExW = zeros(1,25);
Scenario = 0;
NonLinCounter = 0;
NoSolution = 1;

Screensize = get(0,'screensize');
vert = Screensize(4)/4;
horiz = Screensize(3)/5;
fig = {};
for k=1:5
set(groot,'defaultfigureposition',[horiz*(k-1) 100 horiz vert])
fig{k} = figure(k);
end
for k=6:10
set(groot,'defaultfigureposition',[horiz*(k-6) 100+vert horiz vert])
fig{k} = figure(k);
end
for k=11:13
set(groot,'defaultfigureposition',[horiz*(k-11) 100+2*vert horiz vert])
fig{k} = figure(k);
end

FlowMax4 = 4.2*10^3; % In m^3/s.
FlowMax1 = 640; % In m^3/s.
FuelMax = 2*10^6; % In hW.
VolMax4 = 5*10^9; % In m^3.
VolMax1 = 9.5*10^9; % In m^3.
VolMinCond4 = 3.47*10^9; % In m^3.
VolMinCond1 = 1.93*10^9; % In m^3.
TMax = 24*3600; % In s.
DMax = 2*10^6; % In hW.

Vini4 = 0.9; % Initial condition between 0.6 and 1.
Vini1 = 0.8; % Initial condition between 0.2 and 1.
Vol4 = @(i) i*VolMax4;
Vol1 = @(i) i*VolMax1;

Tmax = 1;
Tmin = 0;
NT = 24;
dt = (Tmax-Tmin)/NT;
time = Tmin:dt:Tmax;

IT4 = @(t) 2800; % In m^3/s.
IT1 = @(t) 320; % In m^3/s.
% ===================== VOLUMES ===================== (1,0)
Estiramiento = 0; % 10 por ciento.
Vmax4 = (Vini4*VolMax4+IT4(1)*3600*24)/VolMax4*(1+Estiramiento/100);
Vmin4 = (Vini4*VolMax4+(IT4(1)-FlowMax4)*3600*24)/VolMax4*(1-Estiramiento/100);
Vmax1 = (Vini1*VolMax1+IT1(1)*3600*24)/VolMax1*(1+Estiramiento/100);
Vmin1 = (Vini1*VolMax1+(IT1(1)-FlowMax1)*3600*24)/VolMax1*(1-Estiramiento/100);
% Vmax4 = VolMax4;
% Vmin4 = VolMinCond4;
% Vmax1 = VolMax1;
% Vmin1 = VolMinCond1;
% ===================================================
DV4 = Vmax4-Vmin4;
DV1 = Vmax1-Vmin1;
Num = 300;
Num4 = Num;
Num1 = Num;
% Num4 = 150;
% Num1 = 150;
dV4 = DV4/Num4;
dV1 = DV1/Num1;
V4 = Vmin4:dV4:Vmax4;
V1 = Vmin1:dV1:Vmax1;

% D = DMax*(1/2+(1/2)*sin(8*pi*time)); % In kW.
% D = DMax*2/3*time./time; % In kW.
D = D01012017'*1000;

% === WIND ===>>>
Gk = 1383000; % kW.
theta_W0 = 0.2;
theta_W1 = 0.03;
alpha_W = 0.14;
theta_W = @(t) theta_W0*exp(-theta_W1*t);
pW = WP01012017;
Y(:,1) = Wind_For-SigW;
Y(:,2) = Wind_For;
Y(:,3) = Wind_For+SigW;
Y_hat = [-1,0,1];
DW(:,1) = max(0,D(:)-Y(:,1)*Gk); % Few wind.
DW(:,2) = max(0,D(:)-Y(:,2)*Gk);
DW(:,3) = max(0,D(:)-Y(:,3)*Gk); % High wind.
cw{length(time),length(V1),length(V4)} = {};
% === WIND ===>>>

% === PLOT DEMAND ===>>>
set(0,'CurrentFigure',6);
plot(time,DW(:,1));
hold on
plot(time,DW(:,2));
plot(time,DW(:,3));
title('Demand');
xlabel('Time');
ylabel('hW');
grid on;
pause(0.001);
% === PLOT DEMAND ===>>>

u = zeros(length(time),length(V1),length(V4),3); % 3 wind lines.
Ph4 = zeros(length(time),length(V1),length(V4));
Ph1 = zeros(length(time),length(V1),length(V4));
Pf = zeros(length(time),length(V1),length(V4));
PhiT4 = zeros(length(time),length(V1),length(V4),3);
PhiT1 = zeros(length(time),length(V1),length(V4),3);
PhiF = zeros(length(time),length(V1),length(V4),3);
Ph4_Plot = zeros(length(time),length(V1),length(V4),3);
Ph1_Plot = zeros(length(time),length(V1),length(V4),3);
Pf_Plot = zeros(length(time),length(V1),length(V4),3);
HAM_Plot = zeros(length(time),length(V1),length(V4),3);
uv4 = zeros(length(V1),length(V4));
uv1 = zeros(length(V1),length(V4));

eta = 8.72;
h04 = 5.1;
h01 = 54; % This later will be the level of Baygorria.


% Khe4 = 0.3; % In US$/MWh.
Khe4 = 0.3; % In US$/MWh.
Khe1 = 0.3; % In US$/MWh.
Kf = 130; % In US$/MWh.

H4 = @(V) (-19.8)*(V).^2 + (51.5)*(V) + 3.79;
H1 = @(V) (-3.74)*(V).^2 + (16.7)*(V) + 67.7;

cd4 = 5.3;
d4 = @(tur) cd4*tur;

cd1 = 1.0;
d1 = @(tur) cd1*tur;

PowH1 = @(phihatH1,V1) eta*FlowMax1*phihatH1*(H1(V1)-d1(phihatH1)-h01);
PowH4 = @(phihatH4,V4) eta*FlowMax4*phihatH4*(H4(V4)-d4(phihatH4)-h04);
PowF = @(phihatF) FuelMax*phihatF;

Kh4 = Khe4*eta*(H4(Vini4)-d4(1)-h04)/(3.6*10^6); % Water's value of Salto Grande.
Kh1 = Khe1*eta*(H1(Vini1)-d1(1)-h01)/(3.6*10^6); % Water's value of Bonete.

hatKf = Kf*FuelMax/(3.6*10^6);

hatKh4 = @(uv) FlowMax4*(Kh4-uv/VolMax4);
hatKh1 = @(uv) FlowMax1*(Kh1-uv/VolMax1);

K1_4 = @(V) eta*FlowMax4*(H4(V)-h04);
K1_1 = @(V) eta*FlowMax1*(H1(V)-h01);

K2_4 = -eta*cd4*FlowMax4;
K2_1 = -eta*cd1*FlowMax1;

olK2_4 = -K2_4;
olK2_1 = -K2_1;

Ph_star4 = @(V,uv) (1/(4*olK2_4))*(K1_4(V).^2-(hatKh4(uv)*FuelMax/hatKf).^2);
Ph_star1 = @(V,uv) (1/(4*olK2_1))*(K1_1(V).^2-(hatKh4(uv)*FuelMax/hatKf).^2);

olPh4 = @(V) eta*FlowMax4*(H4(V)-d4(1)-h04);
olPh1 = @(V) eta*FlowMax1*(H1(V)-d1(1)-h01);

phit4 = @(V,Ph) (K1_4(V)-sqrt(K1_4(V)^2-4*olK2_4*Ph))/(2*olK2_4);
phit1 = @(V,Ph) (K1_1(V)-sqrt(K1_1(V)^2-4*olK2_1*Ph))/(2*olK2_1);

CinZ1 = @(V,uv) hatKh1(uv)/K1_1(V) - hatKf/FuelMax;
CinZ4 = @(V,uv) hatKh4(uv)/K1_4(V) - hatKf/FuelMax;

CostPh4 = @(uv1,uv4,V1,V4,D,Ph4) hatKh4(uv4)*(K1_4(V4)-sqrt(K1_4(V4)^2-4*olK2_4*Ph4))/(2*olK2_4) + hatKh1(uv1)*(K1_1(V1)-sqrt(K1_1(V1)^2-4*olK2_1*(D-Ph4)))/(2*olK2_1);
dCostPh4 = @(uv1,uv4,V1,V4,D,Ph4) hatKh4(uv4)./sqrt(K1_4(V4)^2-4*olK2_4.*Ph4) - hatKh1(uv1)./sqrt(K1_1(V1)^2-4*olK2_1.*(D-Ph4));
% CostPh4inZ = @(t,uv1,uv4,V1,V4) -1;
Ph_StarStar = @(D,uv1,uv4,V1,V4) (4*olK2_1*hatKh4(uv4)^2*D + K1_4(V4)^2*hatKh1(uv1)^2 - K1_1(V1)^2*hatKh4(uv4)^2) / (4*(olK2_4*hatKh1(uv1)^2 + olK2_1*hatKh4(uv4)^2));
olP_Star = @(V4) K1_4(V4)^2/(4*olK2_4);
ulP_Star = @(V1,D) (K1_1(V1)^2-4*olK2_1*D)/(-4*olK2_1);

HamToMin = @(uv4,uv1,phit4,phit1,phif,t,cw,Y_hat) TMax*(hatKh4(uv4)*phit4 + hatKh1(uv1)*phit1 + hatKf*phif + IT4(t)*uv4/VolMax4 + IT1(t)*uv1/VolMax1) + (theta_W(t)*(pW(t)-Y_hat*SigW(t))/SigW(t)-Y_hat*d_SigW(t)/SigW(t))*(2*cw(1)*Y_hat+cw(2)) + (theta_W(t)*alpha_W*Y_hat*(1-Y_hat*SigW(t))/SigW(t))*2*cw(1);

A1 = @(V1)eta*FlowMax1*(H1(V1)-h01);
A4 = @(V4)eta*FlowMax4*(H4(V4)-h04);
B1 = eta*cd1*FlowMax1;
B4 = eta*cd4*FlowMax4;
C1 = @(uv1,phi1op,V1) hatKh1(uv1)/(A1(V1)-2*B1*phi1op);
C4 = @(uv4,phi4op,V4) hatKh4(uv4)/(A4(V4)-2*B4*phi4op);
CF = hatKf/FuelMax;

Num_Reduc = 0; % Number of partitions that I reduce in each time step.
Reduction = -Num_Reduc;

for t=length(time):-1:1
    
    disp(['Complete: ',num2str((length(time)-t)*100/length(time)),'%.']);
    uv4 = zeros(length(V1),length(V4));
    uv1 = zeros(length(V1),length(V4));
    Reduction = Reduction + Num_Reduc;
    newStart = 1+Reduction;

    for w=1:3
        
        disp(['Wind Case: ',num2str(w),' out of 3']);

        for i=1+Reduction:length(V1)-Reduction
            for j=1+Reduction:length(V4)-Reduction

                if (i~=1+Reduction) && (j~=1+Reduction) && (i~=length(V1)-Reduction) && (j~=length(V4)-Reduction) && (t~=length(time))
                    uv1(i,j) = (u(t+1,i+1,j,w)-u(t+1,i-1,j,w))/(2*dV1);
                    uv4(i,j) = (u(t+1,i,j+1,w)-u(t+1,i,j-1,w))/(2*dV4);
                else

                    if i==1+Reduction && (t~=length(time))
                        uv1(i,j) = (u(t+1,i+1,j,w)-u(t+1,i,j,w))/(dV1);
                    end
                    if i==length(V1)-Reduction && (t~=length(time))
                        uv1(i,j) = (u(t+1,i,j,w)-u(t+1,i-1,j,w))/(dV1);
                    end
                    if j==1+Reduction && (t~=length(time))
                        uv4(i,j) = (u(t+1,i,j+1,w)-u(t+1,i,j,w))/(dV4);
                    end
                    if j==length(V4)-Reduction && (t~=length(time))
                        uv4(i,j) = (u(t+1,i,j,w)-u(t+1,i,j-1,w))/(dV4);
                    end
                    if (j==1+Reduction || j==length(V4)-Reduction) && (i~=1+Reduction) && (i~=length(V1)-Reduction) && (t~=length(time))
                        uv1(i,j) = (u(t+1,i+1,j,w)-u(t+1,i-1,j,w))/(2*dV1);
                    end
                    if (i==1+Reduction || i==length(V1)-Reduction) && j~=1+Reduction && j~=length(V4)-Reduction && (t~=length(time))
                        uv4(i,j) = (u(t+1,i,j+1,w)-u(t+1,i,j-1,w))/(2*dV4);
                    end

                end
                
                if t==length(time)
                    cw{t,i,j} = [0,0,0]; % I define the polynomial with all ceros.
                end

                
                % ============================ HAMILTONIAN ============================ 

                if CinZ1(V1(i),uv1(i,j))<0 && CinZ4(V4(j),uv4(i,j))<0 && Ph_star1(V1(i),uv1(i,j))>0 && Ph_star4(V4(j),uv4(i,j))>0

                    % ============================ Scenario 1 ============================ 
                    if DW(t,w) <= olPh1(V1(i)) % Scenario 1.
                        Scenario = 1;

                        if Ph_star1(V1(i),uv1(i,j))+Ph_star4(V4(j),uv4(i,j))<=DW(t,w)
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Pf(t,i,j) = DW(t) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j))+Ph_star4(V4(j),uv4(i,j))>DW(t,w)
                            if not((dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),0)>0 && dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),DW(t,w))<0))% && Ph_StarStar(t,uv1(i,j),uv4(i,j),V1(i),V4(j))>0
                                if Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))<=DW(t,w) && Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))>0
                                    Ph4(t,i,j) = Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))>DW(t,w)
                                    Ph1(t,i,j) = 0;
                                    Ph4(t,i,j) = DW(t,w);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))<0
                                    Ph1(t,i,j) = DW(t,w);
                                    Ph4(t,i,j) = 0;
                                    Pf(t,i,j) = 0;
                                else
                                    disp('Error in Scenario 1 (C).');
                                    return;
                                end
                            else
                                FF = ulP_Star(V1(i),DW(t,w)):1000:olP_Star(V4(j));
                                figure;
                                plot(FF,dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),FF));
                                hold on;
                                plot(FF,CostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),FF));
                                vline(0,'r',' ');
                                vline(DW(t,w),'r',' ');
                                vline(Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j)),'b',' ');
                                pause(1.5);   
                                disp('Error in Scenario 1 (B).');
                                return;
                            end
                        else
                            disp('Error in Scenario 1 (A).');
                            return;
                        end

                    end
                    % ============================ Scenario 2 ============================ 
                    if DW(t,w)>olPh1(V1(i)) && DW(t,w)<=olPh4(V4(j))
                        Scenario = 2;

                        if Ph_star1(V1(i),uv1(i,j))+Ph_star4(V4(j),uv4(i,j))<=DW(t,w) && Ph_star1(V1(i),uv1(i,j))<=olPh1(V1(i))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star4(V4(j),uv4(i,j))<=DW(t,w)-olPh1(V1(i)) && Ph_star1(V1(i),uv1(i,j))>olPh1(V1(i))
                            Ph1(t,i,j) = olPh1(V1(i));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j))+Ph_star4(V4(j),uv4(i,j))>DW(t,w) && Ph_star4(V4(j),uv4(i,j))>DW(t,w)-olPh1(V1(i))
                            if not((dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),DW(t,w)-olPh1(V1(i)))>0 && dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),DW(t,w))<0))% && Ph_StarStar(t,uv1(i,j),uv4(i,j),V1(i),V4(j))>0
                                if Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))>=DW(t,w)-olPh1(V1(i)) && Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))<=DW(t,w)
                                    Ph4(t,i,j) = Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))>DW(t,w)
                                    Ph4(t,i,j) = DW(t,w);
                                    Ph1(t,i,j) = 0;
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))<DW(t,w)-olPh1(V1(i))
                                    Ph1(t,i,j) = olPh1(V1(i));
                                    Ph4(t,i,j) = DW(t,w) - Ph1(t,i,j);
                                    Pf(t,i,j) = 0;
                                else
                                    disp('Error in Scenario 2 (C).');
                                    return;
                                end
                            else
                                FF = ulP_Star(V1(i),DW(t,w)):1000:olP_Star(V4(j));
                                figure;
                                plot(FF,dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),FF));
                                hold on;
                                plot(FF,CostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),FF));
                                vline(DW(t,w)-olPh1(V1(i)),'r',' ');
                                vline(DW(t,w),'r',' ');
                                vline(Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j)),'b',' ');
                                pause(1.5);                            
                                disp('Error in Scenario 2 (B).');
                                return;
                            end
                        else
                            disp('Error in Scenario 2 (A).');
                            return;
                        end

                    end
                    % ============================ Scenario 3 ============================ 
                    if DW(t,w)>olPh1(V1(i)) && DW(t,w)>olPh4(V4(j)) && DW(t,w)<=olPh1(V1(i))+olPh4(V4(j))
                        Scenario = 3;

                        if Ph_star1(V1(i),uv1(i,j))+Ph_star4(V4(j),uv4(i,j))<=DW(t,w) && Ph_star4(V4(j),uv4(i,j))<=olPh4(V4(j)) && Ph_star1(V1(j),uv1(i,j))<=olPh1(V1(i))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star4(V4(j),uv4(i,j))<=DW(t,w)-olPh1(V1(i)) && Ph_star1(V1(i),uv1(i,j))>olPh1(V1(i))
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Ph1(t,i,j) = olPh1(V1(i));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star4(V4(j),uv4(i,j))>olPh4(V4(j)) && Ph_star1(V1(i),uv1(i,j))<=DW(t,w)-olPh4(V4(j))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j));
                            Ph4(t,i,j) = olPh4(V4(j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j))+Ph_star4(V4(j),uv4(i,j))>DW(t,w) && Ph_star4(V4(j),uv4(i,j))>DW(t,w)-olPh1(V1(i)) && Ph_star1(V1(i),uv1(i,j))>DW(t,w)-olPh4(V4(j))
                            if not((dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),DW(t,w)-olPh1(V1(i)))>0 && dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),olPh4(V4(j)))<0))% && Ph_StarStar(t,uv1(i,j),uv4(i,j),V1(i),V4(j))>0
                                if Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))>=DW(t,w)-olPh1(V1(i)) && Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))<=olPh4(V4(j))
                                    Ph4(t,i,j) = Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))>olPh4(V4(j))
                                    Ph4(t,i,j) = olPh4(V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j))<DW(t,w)-olPh1(V1(i))
                                    Ph1(t,i,j) = olPh1(V1(i));
                                    Ph4(t,i,j) = DW(t,w) - Ph1(t,i,j);
                                    Pf(t,i,j) = 0;
                                else
                                    disp('Error in Scenario 3 (C).');
                                    return;
                                end
                            else
                                FF = ulP_Star(V1(i),DW(t,w)):1000:olP_Star(V4(j));
                                figure;
                                plot(FF,dCostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),FF));
                                hold on;
                                plot(FF,CostPh4(uv1(i,j),uv4(i,j),V1(i),V4(j),DW(t,w),FF));
                                vline(DW(t,w)-olPh1(V1(i)),'r',' ');
                                vline(olPh4(V4(j)),'r',' ');
                                vline(Ph_StarStar(DW(t,w),uv1(i,j),uv4(i,j),V1(i),V4(j)),'b',' ');
                                pause(1.5);                            
                                disp('Error in Scenario 3 (B).');
                                return;
                            end
                        else
                            disp('Error in Scenario 3 (A).');
                            return;
                        end

                    end
                    % ============================ Scenario 4 ============================ 
                    if DW(t,w)>olPh1(V1(i))+olPh4(V4(j))
                        Scenario = 4;

                        if Ph_star1(V1(i),uv1(i,j))<=olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j))<=olPh4(V4(j))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j))>olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j))<=olPh4(V4(j))
                            Ph1(t,i,j) = olPh1(V1(i));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j))<=olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j))>olPh4(V4(j))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j));
                            Ph4(t,i,j) = olPh4(V4(j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j))>olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j))>olPh4(V4(j))
                            Ph1(t,i,j) = olPh1(V1(i));
                            Ph4(t,i,j) = olPh4(V4(j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        else
                            disp('Error in Scenario 4.');
                            return;
                        end

                    end
                    % ============================ END Scenarios ============================ 
                else
                    disp('Error in contidion dC/dPh(0)<0 or Ph*>0 (first condition).');
                    disp(['Time: ',num2str(100*(length(time)-t)/length(time)),'%']);
                    disp(['CinZ1_test_NEG = ',num2str(CinZ1(V1(i),uv1(i,j)))]);
                    disp(['CinZ4_test_NEG = ',num2str(CinZ4(V4(j),uv4(i,j)))]);
                    disp(['Ph_star1_test_POS = ',num2str(Ph_star1(V1(i),uv1(i,j)))]);
                    disp(['Ph_star4_test_POS = ',num2str(Ph_star4(V4(j),uv4(i,j)))]);
                    return;
                end

                % ============================ HAMILTONIAN ============================

                PhiT4(t,i,j,w) = phit4(V4(j),Ph4(t,i,j));
                PhiT1(t,i,j,w) = phit1(V1(i),Ph1(t,i,j));
                PhiF(t,i,j,w) = Pf(t,i,j)/FuelMax;
                
                Ph1_Plot(t,i,j,w) = Ph1(t,i,j);
                Ph4_Plot(t,i,j,w) = Ph4(t,i,j);
                Pf_Plot(t,i,j,w) = Pf(t,i,j);
                
                if t==length(time)
                    cwH = [0,0,0];
                    HAM(i,j) = 0;
                    HAM_Plot(t,i,j,w) = HAM(i,j);
                else
                    cwH = cw{t+1,i,j};
                    HAM(i,j) = HamToMin(uv4(i,j),uv1(i,j),PhiT4(t+1,i,j,w),PhiT1(t+1,i,j,w),PhiF(t+1,i,j,w),t,cwH,Y_hat(w));
                    HAM_Plot(t,i,j,w) = HAM(i,j);
                end

                % ============================ FMINCON ============================
%                 if NonLinCounter == 1000
% 
%                     while NoSolution ~= 0
% 
%                         Demand = DW(t,w);
%                         VG1 = V1(i);
%                         VG4 = V4(j);
%                         
%                         fun = @(x) HamToMin(uv4(i,j),uv1(i,j),x(2),x(1),x(3),t,cwH,Y_hat(w));
%                         
%                         A = [];
%                         b = [];
%                         Aeq = [];
%                         beq = [];
%                         lb = zeros(1,3);
%                         ub = ones(1,3);
%                         x0 = rand(1,3);
%                         options = optimoptions('fmincon','Display','off');
%                         x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
% 
%                         if abs(x(1)-PhiT1(t,i,j,w))>=0.001 || abs(x(2)-PhiT4(t,i,j,w))>=0.001 || abs(x(3)-PhiF(t,i,j,w))>=0.001
%                             NoSolution = NoSolution + 1;
%                             if NoSolution == 10000
%                                 disp('Error in Optimization');
%                                 disp(['Our solution (sum(P)-DW(t,w))=',num2str(Ph4(t,i,j)+Ph1(t,i,j)+Pf(t,i,j)-DW(t,w))]);
%                                 disp(['MATLAB solution (sum(P)-DW(t,w))=',num2str(PowH1(x(1),V1(i))+PowH4(x(2),V4(j))+PowF(x(3))-DW(t,w))]);
%                                 error2(E2:E2+2) = [abs(x(1)-PhiT1(t,i,j)),abs(x(2)-PhiT4(t,i,j)),abs(x(3)-PhiF(t,i,j))];
%                                 return;
%                             end
%                             else
%                             error2(E2:E2+2) = [abs(x(1)-PhiT1(t,i,j)),abs(x(2)-PhiT4(t,i,j)),abs(x(3)-PhiF(t,i,j))];
%                             E2 = E2+3;
%                             NoSolution = 0;
%                         end
% 
%                     end
%                     NonLinCounter = 0;
%                     NoSolution = 1;
%                 else
%                     NonLinCounter = NonLinCounter + 1;
%                 end
                % ============================ FMINCON ============================

                Matrix_Reduction = ones(length(V1),length(V4));
                if Reduction ~= 0
                    Matrix_Reduction(1:Reduction,:)=NaN;
                    Matrix_Reduction(:,end-Reduction:end)=NaN;
                    Matrix_Reduction(:,1:Reduction)=NaN;
                    Matrix_Reduction(end-Reduction:end,:)=NaN;
                end

                HAM = HAM.*Matrix_Reduction; % We remove the points that we do not care anymore.

                PH1logic = PhiT1(t,i,j) >= 0.99999; % 1=>1.
                PH4logic = PhiT4(t,i,j) >= 0.99999; % 1=>1.
                PHFlogic = 1 - (PhiF(t,i,j) == 0); % 0=>0.

                % === LAMBDA ===>>>
                    if PH1logic==1
                        if PH4logic==1
                            if PHFlogic==0
                                Lambda(t,i,j,w) = CF;
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = CF;
                            end
                        elseif PH4logic~=1
                            if PHFlogic==0
                                Lambda(t,i,j,w) = C4(uv4(i,j),PhiT4(t,i,j),V4(j));
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = C4(uv4(i,j),PhiT4(t,i,j),V4(j));
                                error(E) = C4(uv4(i,j),PhiT4(t,i,j),V4(j)) - CF;
                                E = E+1;
                            end
                        end
                    elseif PH1logic~=1
                        if PH4logic==1
                            if PHFlogic==0
                                Lambda(t,i,j,w) = C1(uv1(i,j),PhiT1(t,i,j),V1(i));
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = C1(uv1(i,j),PhiT1(t,i,j),V1(i));
                                error(E) = C1(uv1(i,j),PhiT1(t,i,j),V1(i)) - CF;
                                E = E+1;
                            end
                        elseif PH4logic~=1
                            if PHFlogic==0
                                Lambda(t,i,j,w) = C1(uv1(i,j),PhiT1(t,i,j),V1(i));
                                error(E) = C1(uv1(i,j),PhiT1(t,i,j),V1(i)) - C4(uv4(i,j),PhiT4(t,i,j),V4(j));
                                E = E+1;
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = CF;
                                error(E) = (2/3)*C1(uv1(i,j),PhiT1(t,i,j),V1(i)) - C4(uv4(i,j),PhiT4(t,i,j),V4(j))/3 - CF/3;
                            end
                        end
                    end
                % === LAMBDA ===>>>

            end
        end
            
        if t~=length(time)
            u(t,:,:,w) = u(t+1,:,:,w) + dt*reshape(HAM,[1,length(V1),length(V4)]); % VERIFICAR.
        end
        
    end
        
                    % === COLORS ===>>>
%                     colormap(parula);
                    colormap(jet);
%                     colormap(hsv);
%                     colormap(hot);
%                     colormap(cool);
%                     colormap(spring);
%                     colormap(summer);
%                     colormap(autumn);
%                     colormap(winter);
%                     colormap(bone);
%                     colormap(copper);
%                     colormap(pink);
%                     colormap(lines);
%                     colormap(colorcube);
%                     colormap(prism);
%                     colormap(flag);
%                     colormap(white);                    
                    % === COLORS ===>>>
                    

    
    
            % === PLOT ===>>>
            set(0,'CurrentFigure',6);
            vline(time(t),'r',' ');
            title(['Demand (time=',num2str(t),')']);

            set(0,'CurrentFigure',1);
            [v1,v4] = meshgrid(V1,V4);
            clf(1);
            hold on;
            for w=1:3
                s0 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(u(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s0.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            zlabel('US$');
            title(['Optimal cost function (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            box;
            view(140,20);
        %     saveas(gcf,'Result_1','epsc');

            set(0,'CurrentFigure',2);
            clf(2);
            hold on;
            for w=1:3
                grad = gradient(squeeze(u(t,newStart:end-Reduction,newStart:end-Reduction,w)));
                s1 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),grad');
                s1.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Gradient of optimal cost function (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',7);
            clf(7);
            hold on;
            for w=1:3
                s2 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(Lambda(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s2.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Lagrangian Value (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',5);
            clf(5);
            hold on;
            for w=1:3
                s3 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(PhiT1(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s3.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Turbined flow 1 (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;
            
            set(0,'CurrentFigure',4);
            clf(4);
            hold on;
            for w=1:3
                s4 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(PhiT4(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s4.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Turbined flow 4 (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',3);
            clf(3);
            hold on;
            for w=1:3
                s5 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(PhiF(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s5.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Fuel control (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',8);
            clf(8);
            hold on;
            for w=1:3
                s6 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(Pf_Plot(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s6.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Fuel Power (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',10);
            clf(10);
            hold on;
            for w=1:3
                s7 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(Ph1_Plot(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s7.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Hydraulic Power 1 (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;
            
            set(0,'CurrentFigure',9);
            clf(9);
            hold on;
            for w=1:3
                s8 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(Ph4_Plot(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s8.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['Hydraulic Power 4 (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;
            
            set(0,'CurrentFigure',11);
            clf(11);
            for w=1:1
                s9 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),uv1(newStart:end-Reduction,newStart:end-Reduction)');
                s9.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['uv1 (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',12);
            clf(12);
            for w=1:1
                s10 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),uv4(newStart:end-Reduction,newStart:end-Reduction)');
                s10.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['uv4 (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;

            set(0,'CurrentFigure',13);
            clf(13);
            hold on;
            for w=1:3
                s11 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(HAM_Plot(t,newStart:end-Reduction,newStart:end-Reduction,w))');
                s11.EdgeColor = 'none';
            end
            xlabel('V1');
            ylabel('V4');
            title(['HAMILTONIAN (time=',num2str(t),')']);
            xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
            ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
            zlim auto
            view(140,20);
            box;
%             pause(1);
            % === PLOT ===>>>
    
    % === CW ===>>>
    if t~=length(time)
        for i=1:length(V1)
            for j=1:length(V4)
                cw{t,i,j} = polyfit(Y_hat,squeeze(u(t,i,j,:))',2);
            end
        end
    end
    % === CW ===>>>
    
    set(groot,'defaultFigureVisible','on') % Show plots ON.
    if t==length(time)
        set(groot,'defaultfigureposition',[0 0 Screensize(3) Screensize(4)]);
        fig{14} = figure(14);
    else
        clf(14);
    end
    for i=1:12
        set(0,'CurrentFigure',i);
        ax{i} = gca;
        set(0,'CurrentFigure',14);
        s{i} = subplot(3,4,i);
        switch i
            case 1
                xlabel('V1');
                ylabel('V4');
                zlabel('US$');
                title(['Optimal cost function (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                box;
                view(140,20);
            case 2
                xlabel('V1');
                ylabel('V4');
                title(['Gradient of optimal cost function (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 3
                xlabel('V1');
                ylabel('V4');
                title(['Fuel control (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 4
                xlabel('V1');
                ylabel('V4');
                title(['Turbined flow 4 (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 5
                xlabel('V1');
                ylabel('V4');
                title(['Turbined flow 1 (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 6
                title(['Demand (time=',num2str(t),')']);
                vline(time(t),'r',' ');
            case 7
                xlabel('V1');
                ylabel('V4');
                title(['Lagrangian Value (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 8
                xlabel('V1');
                ylabel('V4');
                title(['Fuel Power (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 9
                xlabel('V1');
                ylabel('V4');
                title(['Hydraulic Power 4 (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 10
                xlabel('V1');
                ylabel('V4');
                title(['Hydraulic Power 1 (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 11
                xlabel('V1');
                ylabel('V4');
                title(['uv1 (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 12
                xlabel('V1');
                ylabel('V4');
                title(['uv4 (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;  
        end
        figs{i} = get(ax{i},'children');
        copyobj(figs{i},s{i});
    end
    pause(1);
    set(groot,'defaultFigureVisible','off') % Show plots OFF.
    saveas(gcf,['Result_time_',num2str(t)],'epsc');

end