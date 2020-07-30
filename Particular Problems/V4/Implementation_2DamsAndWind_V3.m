close all
clear all
clc

global Demand VG1 VG4

SaveFig = 0;
SaveU = 1;
ShowFig = 0;
% format long;
if ShowFig == 1
    set(groot,'defaultFigureVisible','on');
elseif ShowFig == 0
    set(groot,'defaultFigureVisible','off');
else
    disp('Choose ShowFig between 0 or 1');
end

E = 1; % To evaluate possible errors.
E2 = 1; % To evaluate possible errors in the optimization.
load('D01012017.mat');
load('Wind_Data.mat');
GM = 1383; % MW.
ExW = Wind_Data(1,:)'+0.04;
d_ExW = Wind_Data(2,:)';
SigW = Wind_Data(3,:)'/2;
d_SigW = Wind_Data(4,:)'/2;
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
for k=11:15
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

IT4 = @(t) FlowMax4/2; % In m^3/s.
IT1 = @(t) FlowMax1/2; % In m^3/s.
% ===================== VOLUMES ===================== (1,0)

% To evaluate the central point:
Vmax4 = (Vini4*VolMax4+IT4(1)*3600*24)/VolMax4;
Vmin4 = (Vini4*VolMax4-IT4(1)*3600*24)/VolMax4;
Vmax1 = (Vini1*VolMax1+IT1(1)*3600*24)/VolMax1;
Vmin1 = (Vini1*VolMax1-IT1(1)*3600*24)/VolMax1;

% To save time:
Estiramiento = 0; % 10 por ciento.
% Vmax4 = (Vini4*VolMax4+IT4(1)*3600*24)/VolMax4*(1+Estiramiento/100);
% Vmin4 = (Vini4*VolMax4+(IT4(1)-FlowMax4)*3600*24)/VolMax4*(1-Estiramiento/100);
% Vmax1 = (Vini1*VolMax1+IT1(1)*3600*24)/VolMax1*(1+Estiramiento/100);
% Vmin1 = (Vini1*VolMax1+(IT1(1)-FlowMax1)*3600*24)/VolMax1*(1-Estiramiento/100);

% To see all the possible volumes:
% Vmax4 = VolMax4;
% Vmin4 = VolMinCond4;
% Vmax1 = VolMax1;
% Vmin1 = VolMinCond1;
% ===================================================
DV4 = Vmax4-Vmin4;
DV1 = Vmax1-Vmin1;
Num = 6*2^3;
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
D = D01012017'*1000; % In kW.

% === WIND ===>>>
Gk = 1383000; % kW.
theta_W0 = 0.2;
theta_W1 = 0.03;
alpha_W = 0.14;
theta_W = @(t) theta_W0*exp(-theta_W1*t);
pW = ExW;
alpha = 3.5;
coss = @(y) -alpha*cos(y*pi/6);
Y_hat = coss(0:1:6);
for i=1:length(Y_hat)
    for j=1:length(SigW)
        Y(j,i) = pW(j)+Y_hat(i)*SigW(j)/2;
        DW(j,i) = max(0,D(j)-Y(j,i)*Gk);
    end
end
cw{length(time),length(V1),length(V4)} = {};
% === WIND ===>>>

% === PLOT DEMAND ===>>>
set(0,'CurrentFigure',6);
hold on
for i=1:length(Y_hat)
    plot(time,DW(:,i));
end
title('Effective Demand');
xlabel('Time');
ylabel('kW');
grid on;
pause(0.001);
% === PLOT DEMAND ===>>>

u = zeros(length(time),length(V1),length(V4),length(Y_hat)); % 3 wind lines.
Ph4 = zeros(length(time),length(V1),length(V4));
Ph1 = zeros(length(time),length(V1),length(V4));
Pf = zeros(length(time),length(V1),length(V4));
PhiT4 = zeros(length(time),length(V1),length(V4),length(Y_hat));
PhiT1 = zeros(length(time),length(V1),length(V4),length(Y_hat));
PhiF = zeros(length(time),length(V1),length(V4),length(Y_hat));
Ph4_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
Ph1_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
Pf_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
HAM_Plot = zeros(length(time),length(V1),length(V4),length(Y_hat));
HAM_Plot_HF = zeros(length(time),length(V1),length(V4),length(Y_hat));
HAM_Plot_W = zeros(length(time),length(V1),length(V4),length(Y_hat));
uv4 = zeros(length(V1),length(V4),length(Y_hat));
uv1 = zeros(length(V1),length(V4),length(Y_hat));

eta = 8.72;
h04 = 5.1;
h01 = 54; % This later will be the level of Baygorria.

Khe4 = 1.0; % In US$/MWh.
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
Ph_StarStar = @(D,uv1,uv4,V1,V4) (4*olK2_1*hatKh4(uv4)^2*D + K1_4(V4)^2*hatKh1(uv1)^2 - K1_1(V1)^2*hatKh4(uv4)^2) / (4*(olK2_4*hatKh1(uv1)^2 + olK2_1*hatKh4(uv4)^2));
olP_Star = @(V4) K1_4(V4)^2/(4*olK2_4);
ulP_Star = @(V1,D) (K1_1(V1)^2-4*olK2_1*D)/(-4*olK2_1);

HamToMin = @(uv4,uv1,phit4,phit1,phif,t,cw,Y_hat) TMax*(hatKh4(uv4)*phit4 + hatKh1(uv1)*phit1 + hatKf*phif + IT4(t)*uv4/VolMax4 + IT1(t)*uv1/VolMax1) + (-Y_hat*(theta_W(t)+d_SigW(t)/SigW(t)))*(6*cw(1)*Y_hat^5+5*cw(2)*Y_hat^4+4*cw(3)*Y_hat^3+3*cw(4)*Y_hat^2+2*cw(5)*Y_hat+cw(6)) + (4*theta_W(t)*alpha_W/(SigW(t)^2))*(Y_hat*SigW(t)/2+pW(t))*(1-Y_hat*SigW(t)/2-pW(t))*(6*5*cw(1)*Y_hat^4+5*4*cw(2)*Y_hat^3+4*3*cw(3)*Y_hat^2+3*2*cw(4)*Y_hat+2*cw(5));
HamToMin_HF = @(uv4,uv1,phit4,phit1,phif,t,cw,Y_hat) TMax*(hatKh4(uv4)*phit4 + hatKh1(uv1)*phit1 + hatKf*phif + IT4(t)*uv4/VolMax4 + IT1(t)*uv1/VolMax1);
HamToMin_W = @(uv4,uv1,phit4,phit1,phif,t,cw,Y_hat) (-Y_hat*(theta_W(t)+d_SigW(t)/SigW(t)))*(6*cw(1)*Y_hat^5+5*cw(2)*Y_hat^4+4*cw(3)*Y_hat^3+3*cw(4)*Y_hat^2+2*cw(5)*Y_hat+cw(6)) + (4*theta_W(t)*alpha_W/(SigW(t)^2))*(Y_hat*SigW(t)/2+pW(t))*(1-Y_hat*SigW(t)/2-pW(t))*(6*5*cw(1)*Y_hat^4+5*4*cw(2)*Y_hat^3+4*3*cw(3)*Y_hat^2+3*2*cw(4)*Y_hat+2*cw(5));

A1 = @(V1)eta*FlowMax1*(H1(V1)-h01);
A4 = @(V4)eta*FlowMax4*(H4(V4)-h04);
B1 = eta*cd1*FlowMax1;
B4 = eta*cd4*FlowMax4;
C1 = @(uv1,phi1op,V1) hatKh1(uv1)/(A1(V1)-2*B1*phi1op);
C4 = @(uv4,phi4op,V4) hatKh4(uv4)/(A4(V4)-2*B4*phi4op);
CF = hatKf/FuelMax;

%%% ====
c_op = hatKf/FuelMax;
d1_op = @(uv1) hatKh1(uv1)/(2*olK2_1*sqrt(4*olK2_1));
d4_op = @(uv4) hatKh4(uv4)/(2*olK2_4*sqrt(4*olK2_4));
e1_op = @(V1) (K1_1(V1)^2)/(4*olK2_1);
e4_op = @(V4) (K1_4(V4)^2)/(4*olK2_4);
% Always with 2 dams, j=1 and then i=4.
g4_op = @(uv1,uv4,V1,V4)  e4_op(V4)-e1_op(V1)*(d4_op(uv4)/d1_op(uv1))^2;
h4_op = @(uv1,uv4) (d4_op(uv4)/d1_op(uv1))^2;
%%% ====

Num_Reduc = 0; % Number of partitions that I reduce in each time step.
Reduction = -Num_Reduc;

for t=length(time):-1:1
    
    disp(['Complete: ',num2str((length(time)-t)*100/length(time)),'%.']);
    uv4 = zeros(length(V1),length(V4),length(Y_hat));
    uv1 = zeros(length(V1),length(V4),length(Y_hat));
    Reduction = Reduction + Num_Reduc;
    newStart = 1+Reduction;

    for w=1:length(Y_hat)
        
        disp(['Wind Case: ',num2str(w),' out of ',num2str(length(Y_hat))]);

        for i=1+Reduction:length(V1)-Reduction
            for j=1+Reduction:length(V4)-Reduction

                if (i~=1+Reduction) && (j~=1+Reduction) && (i~=length(V1)-Reduction) && (j~=length(V4)-Reduction) && (t~=length(time))
                    uv1(i,j,w) = (u(t+1,i+1,j,w)-u(t+1,i-1,j,w))/(2*dV1);
                    uv4(i,j,w) = (u(t+1,i,j+1,w)-u(t+1,i,j-1,w))/(2*dV4);
                else

                    if i==1+Reduction && (t~=length(time))
                        uv1(i,j,w) = (u(t+1,i+1,j,w)-u(t+1,i,j,w))/(dV1);
                    end
                    if i==length(V1)-Reduction && (t~=length(time))
                        uv1(i,j,w) = (u(t+1,i,j,w)-u(t+1,i-1,j,w))/(dV1);
                    end
                    if j==1+Reduction && (t~=length(time))
                        uv4(i,j,w) = (u(t+1,i,j+1,w)-u(t+1,i,j,w))/(dV4);
                    end
                    if j==length(V4)-Reduction && (t~=length(time))
                        uv4(i,j,w) = (u(t+1,i,j,w)-u(t+1,i,j-1,w))/(dV4);
                    end
                    if (j==1+Reduction || j==length(V4)-Reduction) && (i~=1+Reduction) && (i~=length(V1)-Reduction) && (t~=length(time))
                        uv1(i,j,w) = (u(t+1,i+1,j,w)-u(t+1,i-1,j,w))/(2*dV1);
                    end
                    if (i==1+Reduction || i==length(V1)-Reduction) && j~=1+Reduction && j~=length(V4)-Reduction && (t~=length(time))
                        uv4(i,j,w) = (u(t+1,i,j+1,w)-u(t+1,i,j-1,w))/(2*dV4);
                    end

                end
                
                if t==length(time)
                    cw{t,i,j} = [0,0,0,0,0,0,0];
                end

                % ============================ HAMILTONIAN ============================ 
                if 1==1%CinZ1(V1(i),uv1(i,j,w))<0 && CinZ4(V4(j),uv4(i,j,w))<0 && Ph_star1(V1(i),uv1(i,j,w))>0 && Ph_star4(V4(j),uv4(i,j,w))>0

                    % ============================ Scenario 1 ============================ 
                    if DW(t,w) <= olPh1(V1(i)) % Scenario 1.
                        Scenario = 1;

                        if Ph_star1(V1(i),uv1(i,j,w))+Ph_star4(V4(j),uv4(i,j,w))<=DW(t,w)
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j,w));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Pf(t,i,j) = DW(t) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j,w))+Ph_star4(V4(j),uv4(i,j,w))>DW(t,w)
                            if not((dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),0)>0 && dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),DW(t,w))<0))% && Ph_StarStar(t,uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>0
                                if Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))<=DW(t,w) && Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>0
                                    Ph4(t,i,j) = Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>DW(t,w)
                                    Ph1(t,i,j) = 0;
                                    Ph4(t,i,j) = DW(t,w);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))<0
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
                                plot(FF,dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),FF));
                                hold on;
                                plot(FF,CostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),FF));
                                vline(0,'r',' ');
                                vline(DW(t,w),'r',' ');
                                vline(Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)),'b',' ');
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

                        if Ph_star1(V1(i),uv1(i,j,w))+Ph_star4(V4(j),uv4(i,j,w))<=DW(t,w) && Ph_star1(V1(i),uv1(i,j,w))<=olPh1(V1(i))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j,w));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star4(V4(j),uv4(i,j,w))<=DW(t,w)-olPh1(V1(i)) && Ph_star1(V1(i),uv1(i,j,w))>olPh1(V1(i))
                            Ph1(t,i,j) = olPh1(V1(i));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j,w))+Ph_star4(V4(j),uv4(i,j,w))>DW(t,w) && Ph_star4(V4(j),uv4(i,j,w))>DW(t,w)-olPh1(V1(i))
                            if 1==1%not((dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),DW(t,w)-olPh1(V1(i)))>0 && dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),DW(t,w))<0))% && Ph_StarStar(t,uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>0
                                if Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>=DW(t,w)-olPh1(V1(i)) && Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))<=DW(t,w)
                                    Ph4(t,i,j) = Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>DW(t,w)
                                    Ph4(t,i,j) = DW(t,w);
                                    Ph1(t,i,j) = 0;
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))<DW(t,w)-olPh1(V1(i))
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
                                plot(FF,dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),FF));
                                hold on;
                                plot(FF,CostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),FF));
                                vline(DW(t,w)-olPh1(V1(i)),'r',' ');
                                vline(DW(t,w),'r',' ');
                                vline(Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)),'b',' ');
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

                        if Ph_star1(V1(i),uv1(i,j,w))+Ph_star4(V4(j),uv4(i,j,w))<=DW(t,w) && Ph_star4(V4(j),uv4(i,j,w))<=olPh4(V4(j)) && Ph_star1(V1(j),uv1(i,j,w))<=olPh1(V1(i))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j,w));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star4(V4(j),uv4(i,j,w))<=DW(t,w)-olPh1(V1(i)) && Ph_star1(V1(i),uv1(i,j,w))>olPh1(V1(i))
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Ph1(t,i,j) = olPh1(V1(i));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star4(V4(j),uv4(i,j,w))>olPh4(V4(j)) && Ph_star1(V1(i),uv1(i,j,w))<=DW(t,w)-olPh4(V4(j))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j,w));
                            Ph4(t,i,j) = olPh4(V4(j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif 1==1%Ph_star1(V1(i),uv1(i,j,w))+Ph_star4(V4(j),uv4(i,j,w))>DW(t,w) && Ph_star4(V4(j),uv4(i,j,w))>DW(t,w)-olPh1(V1(i)) && Ph_star1(V1(i),uv1(i,j,w))>DW(t,w)-olPh4(V4(j))
                            if 1==1%not((dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),DW(t,w)-olPh1(V1(i)))>0 && dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),olPh4(V4(j)))<0))% && Ph_StarStar(t,uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>0
                                if Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>=DW(t,w)-olPh1(V1(i)) && Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))<=olPh4(V4(j))
                                    Ph4(t,i,j) = Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))>olPh4(V4(j))
                                    Ph4(t,i,j) = olPh4(V4(j));
                                    Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
                                    Pf(t,i,j) = 0;
                                elseif Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j))<DW(t,w)-olPh1(V1(i))
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
                                plot(FF,dCostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),FF));
                                hold on;
                                plot(FF,CostPh4(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j),DW(t,w),FF));
                                vline(DW(t,w)-olPh1(V1(i)),'r',' ');
                                vline(olPh4(V4(j)),'r',' ');
                                vline(Ph_StarStar(DW(t,w),uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)),'b',' ');
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

                        if Ph_star1(V1(i),uv1(i,j,w))<=olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j,w))<=olPh4(V4(j))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j,w));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j,w))>olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j,w))<=olPh4(V4(j))
                            Ph1(t,i,j) = olPh1(V1(i));
                            Ph4(t,i,j) = Ph_star4(V4(j),uv4(i,j,w));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j,w))<=olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j,w))>olPh4(V4(j))
                            Ph1(t,i,j) = Ph_star1(V1(i),uv1(i,j,w));
                            Ph4(t,i,j) = olPh4(V4(j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        elseif Ph_star1(V1(i),uv1(i,j,w))>olPh1(V1(i)) && Ph_star4(V4(j),uv4(i,j,w))>olPh4(V4(j))
                            Ph1(t,i,j) = olPh1(V1(i));
                            Ph4(t,i,j) = olPh4(V4(j));
                            Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
                        else
                            disp('Error in Scenario 4.');
                            return;
                        end

                    end
                    % ============================ END Scenarios ============================ 
                else %if with condition
                    disp('Error in contidion dC/dPh(0)<0 or Ph*>0 (first condition).');
                    disp(['Time: ',num2str(100*(length(time)-t)/length(time)),'%']);
                    disp(['CinZ1_test_NEG = ',num2str(CinZ1(V1(i),uv1(i,j,w)))]);
                    disp(['CinZ4_test_NEG = ',num2str(CinZ4(V4(j),uv4(i,j,w)))]);
                    disp(['Ph_star1_test_POS = ',num2str(Ph_star1(V1(i),uv1(i,j,w)))]);
                    disp(['Ph_star4_test_POS = ',num2str(Ph_star4(V4(j),uv4(i,j,w)))]);
                    return;
                end
                    % ========== V1 ==========
                    % ========== V2 ==========
%                 if t~=1
%                     % No restrictions:
%                     Done = 0;
% 
%                     Ph1(t,i,j) = e1_op(V1(i))-(d1_op(uv1(i,j,w))/(2*c_op))^2;
%                     Ph4(t,i,j) = e4_op(V4(j))-(d4_op(uv4(i,j,w))/(2*c_op))^2;
%                     Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j);
%                     if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                         Done = 1;
%                         Scenario = 1; % For debugging.
%                     end
% 
%                     % One restriction:
%                     if not(Done)
%                         CostVal = [];
%                         Optimums = {};
%                         Optimums{1} = 0;
% 
%                         % gD_0:
%                         s = (DW(t,w)-g4_op(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)))/(1+h4_op(uv1(i,j,w),uv4(i,j,w)));
%                         Ph4(t,i,j) = g4_op(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)) + h4_op(uv1(i,j,w),uv4(i,j,w))*s;
%                         Ph1(t,i,j) = DW(t,w)-Ph4(t,i,j);
%                         Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % gD_1:
%                         s = ((DW(t,w)-FuelMax)-g4_op(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)))/(1+h4_op(uv1(i,j,w),uv4(i,j,w)));
%                         Ph4(t,i,j) = g4_op(uv1(i,j,w),uv4(i,j,w),V1(i),V4(j)) + h4_op(uv1(i,j,w),uv4(i,j,w))*s;
%                         Ph1(t,i,j) = (DW(t,w)-FuelMax)-Ph4(t,i,j);
%                         Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         %g1_0:
%                         Ph1(t,i,j) = 0;
%                         Ph4(t,i,j) = e4_op(V4(j)) - (d4_op(uv4(i,j,w))/(2*c_op))^2;
%                         Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         %g1_1:
%                         Ph1(t,i,j) = olPh1(V1(i));
%                         Ph4(t,i,j) = e4_op(V4(j)) - (d4_op(uv4(i,j,w))/(2*c_op))^2;
%                         Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         %g4_0:
%                         Ph4(t,i,j) = 0;
%                         Ph1(t,i,j) = e1_op(V1(i)) - (d1_op(uv1(i,j,w))/(2*c_op))^2;
%                         Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         %g4_1:
%                         Ph4(t,i,j) = olPh4(V4(j));
%                         Ph1(t,i,j) = e1_op(V1(i)) - (d1_op(uv1(i,j,w))/(2*c_op))^2;
%                         Pf(t,i,j) = DW(t,w)-Ph1(t,i,j)-Ph4(t,i,j)-Ph1(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                     end
% 
%                     if length(Optimums)>1
%                         Done = 1;
%                         for y=2:length(Optimums)
%                             if t==length(time)
%                                 cwH = [0,0,0,0,0,0,0];
%                             else
%                                 cwH = cw{t+1,i,j};
%                             end
%                             CostVal(y-1) = HamToMin(uv4(i,j,w),uv1(i,j,w),phit4(V4(j),Optimums{y}(2)),phit1(V1(i),Optimums{y}(1)),Optimums{y}(3)/FuelMax,t,cwH,Y_hat(w));
% 
%                         end
%                         ind_op = find(CostVal == min(CostVal));
%                         Ph1(t,i,j) = Optimums{ind_op+1}(1);
%                         Ph4(t,i,j) = Optimums{ind_op+1}(2);
%                         Pf(t,i,j) = Optimums{ind_op+1}(3);
%                         Scenario = 2; % For debugging.
%                     end
% 
%                     % Two restriction:
%                     if not(Done)
% 
%                         % gD_0+g1_1:
%                         Ph1(t,i,j) = olPh1(V1(i));
%                         Ph4(t,i,j) = DW(t,w) - Ph1(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % gD_0+g4_1:
%                         Ph4(t,i,j) = olPh4(V4(j));
%                         Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % g1_0+gD_1:
%                         Ph1(t,i,j) = 0;
%                         Ph4(t,i,j) = (DW(t,w)-FuelMax) - Ph1(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % g1_0+g4_1:
%                         Ph1(t,i,j) = 0;
%                         Ph4(t,i,j) = olPh4(V4(j));
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % g4_0+gD_1:
%                         Ph4(t,i,j) = 0;
%                         Ph1(t,i,j) = (DW(t,w)-FuelMax) - Ph4(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % g4_0+g1_1:
%                         Ph4(t,i,j) = 0;
%                         Ph1(t,i,j) = olPh1(V1(i));
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % gD_1+g1_1:
%                         Ph1(t,i,j) = olPh1(V1(i));
%                         Ph4(t,i,j) = (DW(t,w)-FuelMax) - Ph1(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % gD_1+g4_1:
%                         Ph4(t,i,j) = olPh4(V4(j));
%                         Ph1(t,i,j) = (DW(t,w)-FuelMax) - Ph4(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % g1_1+g4_1:
%                         Ph4(t,i,j) = olPh4(V4(j));
%                         Ph1(t,i,j) = olPh1(V1(i));
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % gD_0+g1_0:
%                         Ph1(t,i,j) = 0;
%                         Ph4(t,i,j) = DW(t,w) - Ph1(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % gD_0+g4_0:
%                         Ph4(t,i,j) = 0;
%                         Ph1(t,i,j) = DW(t,w) - Ph4(t,i,j);
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                         % g1_0+g4_0:
%                         Ph4(t,i,j) = 0;
%                         Ph1(t,i,j) = 0;
%                         Pf(t,i,j) = DW(t,w) - Ph1(t,i,j) - Ph4(t,i,j);
%                         if (FuelMax>=Pf(t,i,j)) && (Pf(t,i,j)>=0) && (olPh1(V1(i))>=Ph1(t,i,j)) && (olPh4(V4(j))>=Ph4(t,i,j)) && (Ph1(t,i,j)>=0) && (Ph4(t,i,j)>=0)
%                             Optimums{end+1} = [Ph1(t,i,j),Ph4(t,i,j),Pf(t,i,j)];
%                         end
% 
%                     end
% 
%                     if not(Done)
%                         if length(Optimums)>1
%                             for y=2:length(Optimums)
%                                 if t==length(time)
%                                     cwH = [0,0,0,0,0,0,0];
%                                 else
%                                     cwH = cw{t+1,i,j};
%                                 end
%                                 CostVal(y-1) = HamToMin(uv4(i,j,w),uv1(i,j,w),phit4(V4(j),Optimums{y}(2)),phit1(V1(i),Optimums{y}(1)),Optimums{y}(3)/FuelMax,t,cwH,Y_hat(w));
%                             end
%                             ind_op = find(CostVal == min(CostVal));
%                             Ph1(t,i,j) = Optimums{ind_op+1}(1);
%                             Ph4(t,i,j) = Optimums{ind_op+1}(2);
%                             Pf(t,i,j) = Optimums{ind_op+1}(3);
%                             Scenario = 3; % For debugging.
%                         else
%                             disp('Error in our Optimization, it neves finishes');
%                         end
%                     end
%                 end
                % ========== V2 ==========
                % ============================ HAMILTONIAN ============================

                PhiT4(t,i,j,w) = phit4(V4(j),Ph4(t,i,j));
                PhiT1(t,i,j,w) = phit1(V1(i),Ph1(t,i,j));
                PhiF(t,i,j,w) = Pf(t,i,j)/FuelMax;
                
                Ph1_Plot(t,i,j,w) = Ph1(t,i,j);
                Ph4_Plot(t,i,j,w) = Ph4(t,i,j);
                Pf_Plot(t,i,j,w) = Pf(t,i,j);
                
                if t==length(time)
                    cwH = [0,0,0,0,0,0,0];
                    HAM(i,j) = 0;
                    HAM_Plot(t,i,j,w) = HAM(i,j);
                    HAM_Plot_HF(t,i,j,w) = HAM(i,j);
                    HAM_Plot_W(t,i,j,w) = HAM(i,j);
                else
                    cwH = cw{t+1,i,j};
                    HAM(i,j) = HamToMin(uv4(i,j,w),uv1(i,j,w),PhiT4(t+1,i,j,w),PhiT1(t+1,i,j,w),PhiF(t+1,i,j,w),t+1,cwH,Y_hat(w));
                    HAM_Plot(t,i,j,w) = HAM(i,j);
                    HAM_Plot_HF(t,i,j,w) = HamToMin_HF(uv4(i,j,w),uv1(i,j,w),PhiT4(t+1,i,j,w),PhiT1(t+1,i,j,w),PhiF(t+1,i,j,w),t+1,cwH,Y_hat(w));
                    HAM_Plot_W(t,i,j,w) = HamToMin_W(uv4(i,j,w),uv1(i,j,w),PhiT4(t+1,i,j,w),PhiT1(t+1,i,j,w),PhiF(t+1,i,j,w),t+1,cwH,Y_hat(w));
                end

                % ============================ FMINCON ============================
                if t~=1
                    if NonLinCounter == 10000

                        while NoSolution ~= 0

                            Demand = DW(t,w);
                            VG1 = V1(i);
                            VG4 = V4(j);

                            fun = @(x) HamToMin(uv4(i,j,w),uv1(i,j,w),x(2),x(1),x(3),t,cwH,Y_hat(w));

                            A = [];
                            b = [];
                            Aeq = [];
                            beq = [];
                            lb = zeros(1,3);
                            ub = ones(1,3);
                            x0 = rand(1,3);
                            options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',5000,'FunctionTolerance',1.0000e-07);
                            x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

                            if abs(x(1)-PhiT1(t,i,j,w))>=0.001 || abs(x(2)-PhiT4(t,i,j,w))>=0.001 || abs(x(3)-PhiF(t,i,j,w))>=0.001 %|| abs((HamToMin(uv4(i,j,w),uv1(i,j,w),x(2),x(1),x(3),t,cwH,Y_hat(w))-HamToMin(uv4(i,j,w),uv1(i,j,w),PhiT4(t,i,j,w),PhiT1(t,i,j,w),PhiF(t,i,j,w),t,cwH,Y_hat(w)))/HamToMin(uv4(i,j,w),uv1(i,j,w),x(2),x(1),x(3),t,cwH,Y_hat(w)))>=0.001
                                NoSolution = NoSolution + 1;
                                if NoSolution == 100
                                    disp('Error in Optimization');
                                    disp(['Our solution (sum(P)-DW(t,w))=',num2str(Ph4(t,i,j)+Ph1(t,i,j)+Pf(t,i,j)-DW(t,w))]);
                                    disp(['Our solution H = ',num2str(HamToMin(uv4(i,j,w),uv1(i,j,w),PhiT4(t,i,j,w),PhiT1(t,i,j,w),PhiF(t,i,j,w),t,cwH,Y_hat(w)))]);
                                    disp(['MATLAB solution (sum(P)-DW(t,w))=',num2str(PowH1(x(1),V1(i))+PowH4(x(2),V4(j))+PowF(x(3))-DW(t,w))]);
                                    disp(['MATLAB solution H = ',num2str(HamToMin(uv4(i,j,w),uv1(i,j,w),x(2),x(1),x(3),t,cwH,Y_hat(w)))]);
                                    error2(E2:E2+2) = [abs(x(1)-PhiT1(t,i,j,w)),abs(x(2)-PhiT4(t,i,j,w)),abs(x(3)-PhiF(t,i,j,w))];
                                    return;
                                end
                            else
                                error2(E2:E2+2) = [abs(x(1)-PhiT1(t,i,j,w)),abs(x(2)-PhiT4(t,i,j,w)),abs(x(3)-PhiF(t,i,j,w))];
                                E2 = E2+3;
                                NoSolution = 0;
                            end
                        end
                        NonLinCounter = 0;
                        NoSolution = 1;
                    else
                        NonLinCounter = NonLinCounter + 1;
                    end
                end
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
                                Lambda(t,i,j,w) = C4(uv4(i,j,w),PhiT4(t,i,j),V4(j));
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = C4(uv4(i,j,w),PhiT4(t,i,j),V4(j));
                                error(E) = C4(uv4(i,j,w),PhiT4(t,i,j),V4(j)) - CF;
                                E = E+1;
                            end
                        end
                    elseif PH1logic~=1
                        if PH4logic==1
                            if PHFlogic==0
                                Lambda(t,i,j,w) = C1(uv1(i,j,w),PhiT1(t,i,j),V1(i));
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = C1(uv1(i,j,w),PhiT1(t,i,j),V1(i));
                                error(E) = C1(uv1(i,j,w),PhiT1(t,i,j),V1(i)) - CF;
                                E = E+1;
                            end
                        elseif PH4logic~=1
                            if PHFlogic==0
                                Lambda(t,i,j,w) = C1(uv1(i,j,w),PhiT1(t,i,j),V1(i));
                                error(E) = C1(uv1(i,j,w),PhiT1(t,i,j),V1(i)) - C4(uv4(i,j,w),PhiT4(t,i,j),V4(j));
                                E = E+1;
                            elseif PHFlogic~=0
                                Lambda(t,i,j,w) = CF;
                                error(E) = (2/3)*C1(uv1(i,j,w),PhiT1(t,i,j),V1(i)) - C4(uv4(i,j,w),PhiT4(t,i,j),V4(j))/3 - CF/3;
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
%     colormap(parula);
%     colormap(jet);
%     colormap(hsv);
%     colormap(hot);
%     colormap(cool);
%     colormap(spring);
%     colormap(summer);
%     colormap(autumn);
%     colormap(winter);
%     colormap(bone);
%     colormap(copper);
%     colormap(pink);
%     colormap(lines);
%     colormap(colorcube);
%     colormap(prism);
%     colormap(flag);
%     colormap(white);                    
    % === COLORS ===>>>

    % === PLOT ===>>>
    hm = 3;

    set(0,'CurrentFigure',6);
    vline(time(t),'r',' ');
    title(['Demand (time=',num2str(t),')']);

    set(0,'CurrentFigure',1);
    [v1,v4] = meshgrid(V1,V4);
    clf(1);
    hold on;
    for w=1:hm:length(Y_hat)
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

    set(0,'CurrentFigure',2);
    clf(2);
    hold on;
    for w=1:hm:length(Y_hat)
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
    for w=1:hm:length(Y_hat)
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
    for w=1:3:length(Y_hat)
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
    for w=1:hm:length(Y_hat)
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
    for w=1:hm:length(Y_hat)
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
    for w=1:hm:length(Y_hat)
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
    for w=1:hm:length(Y_hat)
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
    for w=1:hm:length(Y_hat)
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
    hold on;
    for w=1:hm:length(Y_hat)
        s9 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(uv1(newStart:end-Reduction,newStart:end-Reduction,w))');
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
    hold on;
    for w=1:hm:length(Y_hat)
        s10 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(uv4(newStart:end-Reduction,newStart:end-Reduction,w))');
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
    for w=1:hm:length(Y_hat)
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

    set(0,'CurrentFigure',14);
    clf(14);
    hold on;
    for w=1:hm:length(Y_hat)
        s12 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(HAM_Plot_HF(t,newStart:end-Reduction,newStart:end-Reduction,w))');
        s12.EdgeColor = 'none';
    end
    xlabel('V1');
    ylabel('V4');
    title(['HAMILTONIAN HF (time=',num2str(t),')']);
    xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
    ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
    zlim auto
    view(140,20);
    box;

    set(0,'CurrentFigure',15);
    clf(15);
    hold on;
    for w=1:hm:length(Y_hat)
        s13 = surf(v1(newStart:end-Reduction,newStart:end-Reduction),v4(newStart:end-Reduction,newStart:end-Reduction),squeeze(HAM_Plot_W(t,newStart:end-Reduction,newStart:end-Reduction,w))');
        s13.EdgeColor = 'none';
    end
    xlabel('V1');
    ylabel('V4');
    title(['HAMILTONIAN W (time=',num2str(t),')']);
    xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
    ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
    zlim auto
    view(140,20);
    box;

    pause(1);
    % === PLOT ===>>>
    
    % === CW ===>>>
    if t~=length(time)
        for i=1:length(V1)
            for j=1:length(V4)
                cw{t,i,j} = polyfit(Y_hat,squeeze(u(t,i,j,:))',6);
            end
        end
    end
    % === CW ===>>>
    
    % === PLOT 16 (all together) ===>>>
    if t==length(time)
        set(groot,'defaultfigureposition',[0 0 Screensize(3) Screensize(4)]);
        fig{16} = figure(16);
    else
        clf(16);
    end
    for i=1:15
        set(0,'CurrentFigure',i);
        ax{i} = gca;
        set(0,'CurrentFigure',16);
        ss{i} = subplot(3,5,i);
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
            case 13
                xlabel('V1');
                ylabel('V4');
                title(['HAMILTONIAN (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 14
                xlabel('V1');
                ylabel('V4');
                title(['HAMILTONIAN HF (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
            case 15
                xlabel('V1');
                ylabel('V4');
                title(['HAMILTONIAN W (time=',num2str(t),')']);
                xlim([Vmin1+dV1*Reduction Vmax1-dV1*Reduction]);
                ylim([Vmin4+dV4*Reduction Vmax4-dV4*Reduction]);
                zlim auto
                view(140,20);
                box;
        end
        figs{i} = get(ax{i},'children');
        copyobj(figs{i},ss{i});
    end
    pause(0.5);
    if SaveFig == 1
    saveas(gcf,['Result_time_',num2str(t)],'epsc');
    end
    % === PLOT 16 (all together) ===>>>
    
    if t==length(time)
        set(groot,'defaultfigureposition',[0 0 Screensize(3)/2 Screensize(4)/2]);
        fig{17} = figure(17);
        sV1 = floor(length(V1)/5);
        sV4 = floor(length(V4)/5);
    end
    set(0,'CurrentFigure',17);
    title(['Cost Vs Wind (time=',num2str(t),')']);
    xlabel('Normalized wind (-2,2)');
    ylabel('US$');
    grid on;
    hold on;
    
    plot(Y_hat,squeeze(u(t,sV1*2,sV4*2,:)),'*')
    plot(-alpha:0.1:alpha,polyval(cw{t,sV1*2,sV4*2},-alpha:0.1:alpha));
%     plot(Y_hat,squeeze(u(t,sV1,sV4,:)),'*')
%     plot([-alpha:0.1:alpha],polyval(cw{t,sV1,sV4},[-alpha:0.1:alpha]));
%     plot(Y_hat,squeeze(u(t,sV1*4,sV4,:)),'*')
%     plot([-alpha:0.1:alpha],polyval(cw{t,sV1*4,sV4},[-alpha:0.1:alpha]));
%     plot(Y_hat,squeeze(u(t,sV1,sV4*4,:)),'*')
%     plot([-alpha:0.1:alpha],polyval(cw{t,sV1,sV4*4},[-alpha:0.1:alpha]));
%     plot(Y_hat,squeeze(u(t,sV1*4,sV4*4,:)),'*')
%     plot([-alpha:0.1:alpha],polyval(cw{t,sV1*4,sV4*4},[-alpha:0.1:alpha]));
%     legend('Low V1 Low V4','Low V1 Low V4','High V1 Low V4','High V1 Low V4','Low V1 High V4','Low V1 High V4','High V1 High V4','High V1 High V4');
    if SaveFig == 1
        saveas(gcf,['Cost_Wind_',num2str(t)],'epsc');
    end
    pause(0.5);
    
    %======== PLOT 18 ========
    if t==length(time)
        set(groot,'defaultfigureposition',[0 0 Screensize(3)/2 Screensize(4)/2]);
        fig{18} = figure(18);
        hold on;
    end

    set(0,'CurrentFigure',18);
    minA = min(min(Y))*0.9;
    maxA = max(max(Y))*1.1;
    
    C = Y_hat*SigW(t)/2+ExW(t);
    D = [minA:(min(C)-minA)/10:min(C),C(2:end-1),max(C):(maxA-max(C))/10:maxA];
    F = 2*(D-ExW(t))/SigW(t);

    for i=1:4
        ss{i} = subplot(2,2,i);
        hold on;
        switch i
            case 1
                xlabel('Wind (0,1)');
                ylabel('US$');
                title('Low V1 Low V4');
                Sol = polyval(cw{t,sV1,sV4},F);
                plot(D(11:17),squeeze(u(t,sV1,sV4,:)),'*');
                plot(D,Sol);
                xlim([minA,maxA]);
                ylim([0,max(squeeze(u(t,sV1,sV4,:)))+1]);
            case 2
                xlabel('Wind (0,1)');
                ylabel('US$');
                title('High V1 Low V4');
                Sol = polyval(cw{t,sV1*4,sV4},F);
                plot(D(11:17),squeeze(u(t,sV1*4,sV4,:)),'*');
                plot(D,Sol);
                xlim([minA,maxA]);
                ylim([0,max(squeeze(u(t,sV1*4,sV4,:)))+1]);
            case 3
                xlabel('Wind (0,1)');
                ylabel('US$');
                title('Low V1 High V4');
                Sol = polyval(cw{t,sV1,sV4*4},F);
                plot(D(11:17),squeeze(u(t,sV1,sV4*4,:)),'*');
                plot(D,Sol);
                xlim([minA,maxA]);
                ylim([0,max(squeeze(u(t,sV1,sV4*4,:)))+1]);
            case 4
                xlabel('Wind (0,1)');
                ylabel('US$');
                title('High V1 High V4');
                Sol = polyval(cw{t,sV1*4,sV4*4},F);
                plot(D(11:17),squeeze(u(t,sV1*4,sV4*4,:)),'*');
                plot(D,Sol);
                xlim([minA,maxA]);
                ylim([0,max(squeeze(u(t,sV1*4,sV4*4,:)))+1]);
        end
    end
    if SaveFig == 1
        saveas(gcf,['All_Winds_',num2str(t)],'epsc');
    end
    pause(0.5)
    %======== PLOT 18 ========

end

if SaveU == 1
    save(['u_',num2str(length(V1)),'x',num2str(length(V4)),'.mat',],'u');
end