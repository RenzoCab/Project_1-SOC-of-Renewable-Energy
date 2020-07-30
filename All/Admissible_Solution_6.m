function [Final_Cost] = Admissible_Solution_6(name,WhatToDo,WarmStart,UseBattery,...
    ComputePlots,NewPlots,SavePlots,Expan,Expan_Time,UseLamb21,UseLamb32,...
    Lambda21,Lambda32,Day,GradNum)

%     Lambda21 = Lambda21/1e6; % We pass from USD/hm^3 to USD/m^3.
%     Lambda32 = Lambda32/1e6; % We pass from USD/hm^3 to USD/m^3.
    % GradNum is in {0,1,...,n-2} where n is Expan_Time.
    Derivatives = [1,1,1,1]; % For debugging.
    
    if 0 == exist(['Simulation_',name],'dir') && WhatToDo ~= 6 && SavePlots == 1
        mkdir(['Simulation_',name]);
    end

    load('MaxFlowCoeff.mat');
    load([pwd '/Historical/Day_',num2str(Day),'_2019.mat']);
    
    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',1e4,...
        'FunctionTolerance',1e-10,'MaxIterations',1e3,'MaxFunEvals',1e7,'Algorithm','sqp',...
        'ConstraintTolerance',1e-5,'StepTolerance',1e-15);

    Plot_11 = 1;
    Plot_21 = 15;
    Plot_12 = 101;
    Plot_22 = 127;
    if NewPlots == 1
        close all;
        if WhatToDo == 1 || WhatToDo == 3
            NeededPlots(Plot_11,Plot_21,0);
        end
        if WhatToDo == 2 || WhatToDo >= 3
            NeededPlots(Plot_12,Plot_22,0);
        end
        figure(11);
        figure(1);
        figure(15);
    elseif NewPlots ~= 0
        disp('Choose NewPlots 0 or 1.');
        return;
    end
    
    MWhTokW = 3.6e6; % From MWh to kW (energy in 1 h to power per second).
    m3Tohm3 = 1e6;  % From m^3 to hm^3.
%     Termicas = {'CTR','Motores Batlle','PTA','PTB','Failure'};
    Termicas = {'Motores Batlle','PTA','PTB','CTR','Failure'};
    Dams = {'Bonete','Baygorria','Palmar','Salto Grande'};
%     MaxTermicas = [Matrix{2}(2),Matrix{2}(4),Matrix{2}(6),Matrix{2}(8),1000]*1000; % kW.
    MaxTermicas = [Matrix{2}(4),Matrix{2}(6),Matrix{2}(8),Matrix{2}(2),1000]*1000; % kW.
    InfTermicas = [0,0,0,0,0];
%     CostTermicas = [Matrix{2}(1),Matrix{2}(3),Matrix{2}(5),Matrix{2}(7),600]; % USD/MWh.
    CostTermicas = [Matrix{2}(3),Matrix{2}(5),Matrix{2}(7),Matrix{2}(1),600]; % USD/MWh.
    CostTermicasNorm = CostTermicas/MWhTokW; % USD/kW.

    NT = 2^Expan_Time; % Discretizations in time.
    NV1 = 2^Expan; % Discretizations of Bonete.
    NV3 = 2^Expan; % Discretizations of Palmar.
    NV4 = 2^Expan; % Discretizations of Salto Grande.
    if UseBattery == 1
        NA = 2^(Expan+3); % Discretizations of Salto Grande.
    else
        NA = 2;
    end
    
    Pc2 = 13;
    Pc3 = 14;

    % ==================== Time ====================>>>

    Tmax = 1; % Simulation final time.
    Tmin = 0; % Simulation initial time.
    dt = (Tmax-Tmin)/NT;
    time = Tmin:dt:Tmax; % Simulation time.

    % ==================== Parameters (Demand, Water and Fuel) ====================>>
       
    H1 = @(V1) (-3.77)*(V1).^2 + (15.7)*(V1) + 69.6; % Height of the water as a function of the normal volume (Bonete).
    H2 = 54; % We fix the level of Baygorria to its nominal level.
    H3 = @(V3) (-7.50)*(V3).^2 + (23.7)*(V3) + 25.8; % Height of the water as a function of the normal volume (Palmar).
    H4 = @(V4) (-21.2)*(V4).^2 + (54.8)*(V4) + 1.92; % Height of the water (m) as a function of the normal volume (Salto Grande).
    VolSG = @(H4) 0.0031*H4.^2 - 0.15*H4 + 2.4;

    h04_Vec = Matrix{8}(:,2);
    h03_Vec = Matrix{8}(:,1);
    h04_Vec = interp1(0:1/(length(h04_Vec)-1):1,h04_Vec,time);
    h03_Vec = interp1(0:1/(length(h03_Vec)-1):1,h03_Vec,time);
    h04 = @(t) h04_Vec(t);
    h03 = @(t) h03_Vec(t);
    
    dBB = Matrix{8}(1,7);
    dBP = Matrix{8}(1,8);
    
    dH1 = @(V1) H1(V1) - H2 - dBB;
    dH2 = @(V3) H2 - H3(V3) - dBP;
    dH3 = @(V3,t) H3(V3) - h03(t);
    dH4 = @(V4,t) H4(V4) - h04(t);
    
    Co1 = MaxFlowCoeff{1};
    Co2 = MaxFlowCoeff{2};
    Co3 = MaxFlowCoeff{3};

    FlowMax1 = @(V1) Co1(1)*dH1(V1).^2 + Co1(2)*dH1(V1) + Co1(3); % In m^3/s.
    FlowMax2 = @(V3) (dH2(V3) <= Co2{3}).*(Co2{1}(1)*dH2(V3) + Co2{1}(2)) + ...
        (dH2(V3) > Co2{3}).*(Co2{2}(1)*dH2(V3) + Co2{2}(2)); % In m^3/s.
    FlowMax3 = @(V3,t) (dH3(V3,t) <= Co3{3}).*(Co3{1}(1)*dH3(V3,t) + Co3{1}(2)) + ...
        (dH3(V3,t) > Co3{3}).*(Co3{2}(1).*dH3(V3,t).^2 + Co3{2}(2).*dH3(V3,t) + Co3{2}(3)); % In m^3/s.
    FlowMax4 = 4410; % In m^3/s.
    
    Vols = [0:0.01:1];
    RealMaxFlow1 = max(FlowMax1(Vols)); % In m^3/s.
    RealMaxFlow2 = max(FlowMax2(Vols)); % In m^3/s.
    RealMaxFlow3 = max(FlowMax3(Vols,find(h03_Vec==max(h03_Vec)))); % In m^3/s.
    RealMaxFlow4 = FlowMax4; % In m^3/s.
    
    MaxTot1 = RealMaxFlow1*2;
    MaxTot2 = RealMaxFlow2*10;
    MaxTot3 = RealMaxFlow3*2;
    MaxTot4 = RealMaxFlow4*2;
    
    MaxSpill1 = @(V1) MaxTot1 - FlowMax1(V1);
    MaxSpill2 = @(V3) MaxTot2 - FlowMax2(V3);
    MaxSpill3 = @(V3,t) MaxTot3 - FlowMax3(V3,t);
    MaxSpill4 = MaxTot4 - RealMaxFlow4;

    MaxV2 = MaxTot1;
    MaxV3 = MaxTot2;

    VolMax1 = 10.7e9; % In m^3.
    VolMax2 = 1; % We do not care about this while H2 is a constant.
    VolMax3 = 3.53e9; % In m^3.
    VolMax4 = 5.18e9; % In m^3.

    VolMin1 = 1.85e9; % In m^3. From tables at 70 m.
    VolMin3 = 1.36e9; % In m^3. From tables at 34 m.
    VolMin4 = 3.65e9; % In m^3. From SimSEE at 30 m.

    TMax = 24*3600; % In s. Number of seconds in a day.
    
    Inflow1 = Matrix{4}(2)*1e6/TMax; % Flow in m^3/s.
    Inflow2 = Matrix{5}(2)*1e6/TMax; % Flow in m^3/s.
    Inflow3 = Matrix{6}(2)*1e6/TMax; % Flow in m^3/s.
    Inflow4 = Matrix{7}(1); % Flow in m^3/s.
    IT1 = @(t) Inflow1; % In m^3/s.
    IT2 = @(t) Inflow2; % In m^3/s.
    IT3 = @(t) Inflow3; % In m^3/s.
    IT4 = @(t) Inflow4/2; % In m^3/s.

    f1 = @(t,V1,T,S) (TMax/VolMax1)*(IT1(t)-FlowMax1(V1)*T-MaxSpill1(V1)*S);
    f3 = @(t,V3,T,S,V) (TMax/VolMax3)*(IT3(t)-FlowMax3(V3,t)*T-MaxSpill3(V3,t)*S+MaxV3*V);
    f4 = @(t,V4,T,S) (TMax/VolMax4)*(IT4(t)-FlowMax4*T-MaxSpill4*S);

    DMax = 2.5e6; % In kW. The maximum possible value of the demand.

    Vini1 = Matrix{4}(1)/VolMax1*1000000;
    Vini3 = Matrix{6}(1)/VolMax3*1000000;
    Vini4 = VolSG(Matrix{7}(2));

    t21 = 1/4; % Delay between Bonete and Baygorria.
    Z21 = floor(length(time)*t21); % Discrete delay between Bonete and Baygorria.
    t32 = 10/24; % Delay between Baygorria and Palmar.
    Z32 = floor(length(time)*t32); % Discrete delay between Baygorria and Palmar.
    
    Phi_V2 = Matrix{9}(:,1)/MaxV2;
    Phi_V3 = Matrix{9}(:,2)/MaxV3;
    Phi_V2 = interp1(linspace(0,1,length(Phi_V2)),Phi_V2,linspace(0,1,Z21));
    Phi_V3 = interp1(linspace(0,1,length(Phi_V3)),Phi_V3,linspace(0,1,Z32));

    if UseLamb21 == 0 || (UseLamb21 == 1 && length(Lambda21) == 1)
        Lambda21 = zeros(1,length(time)+Z21);
    end
    if UseLamb32 == 0 || (UseLamb32 == 1 && length(Lambda32) == 1)
        Lambda32 = zeros(1,length(time)+Z32);
    end
    Lambda21 = interp1(0:1/(length(Lambda21)-1):1,Lambda21,0:1/(length(time)+Z21-1):1);
    Lambda32 = interp1(0:1/(length(Lambda32)-1):1,Lambda32,0:1/(length(time)+Z32-1):1);
    if ComputePlots == 1
        set(0,'CurrentFigure',11); clf(11); set(gcf,'name','Lagrangian Multipliers');
        lamTime = linspace(0,1+t32,length([Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]));
        plot(lamTime,[Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]*1e6);
        hold on; grid on;
        plot(lamTime,[Lambda32,zeros(1,length(Lambda21)-length(Lambda32))]*1e6);
        title('Approximated Lagrangian Multipliers');
        ylabel('USD/hm$^3$','Interpreter','latex');
        xlabel('Time');
        Leg = legend('$\hat{\lambda}_{2,1}(t)$','$\hat{\lambda}_{3,2}(t)$');
        set(Leg, 'interpreter', 'latex');
        xlim([0 max(lamTime)]);
    elseif ComputePlots ~= 0
        disp('Choose ComputePlots between 0 or 1.');
        return;
    end

    NormCost = 1e5; % The normalization is over 100.000 dollars.

    Vmax1 = min((Vini1*VolMax1+(MaxTot1)*TMax),VolMax1)/VolMax1;
    Vmin1 = max((Vini1*VolMax1-(MaxTot1)*TMax),VolMin1)/VolMax1;
    Vmax3 = min((Vini3*VolMax3+(RealMaxFlow3)*TMax),VolMax3)/VolMax3;
    Vmin3 = max((Vini3*VolMax3-(RealMaxFlow3)*TMax),VolMin3)/VolMax3;
    Vmax4 = min((Vini4*VolMax4+(RealMaxFlow4)*TMax),VolMax4)/VolMax4;
    Vmin4 = max((Vini4*VolMax4-(RealMaxFlow4)*TMax),VolMin4)/VolMax4;
    
    dV1 = (Vmax1-Vmin1)/NV1;
    dV3 = (Vmax3-Vmin3)/NV3;
    dV4 = (Vmax4-Vmin4)/NV4;

    V1 = Vmin1:dV1:Vmax1; % Discretized volume of Bonete.
    V3 = Vmin3:dV3:Vmax3; % Discretized volume of Palmar.
    V4 = Vmin4:dV4:Vmax4; % Discretized volume of Salto Grande.

    [vplot1_4,vplot3_4] = meshgrid(V1,V3); % For plotting surfaces.
    [vplot1_3,vplot4_3] = meshgrid(V1,V4); % For plotting surfaces.
    [vplot3_1,vplot4_1] = meshgrid(V3,V4); % For plotting surfaces.
    
    D = Matrix{1}(:,13)'*1000;%+400e3;
    
    Export = (Matrix{1}(:,14)+Matrix{1}(:,15)+Matrix{1}(:,16))*1000;
    WindPower = Matrix{1}(:,6)'*1000;
    SolarPower = Matrix{1}(:,7)'*1000;
    BioMasPower = Matrix{1}(:,9)'*1000;
    timeDemand = 0:1/(length(D)-1):1;
    D = interp1(timeDemand,D,time); %+DMax*0.15*sin(4*pi*time);
    WindPower = interp1(timeDemand,WindPower,time);
    SolarPower = interp1(timeDemand,SolarPower,time);
    BioMasPower = interp1(timeDemand,BioMasPower,time);
    Export = interp1(timeDemand,Export,time);
    

    % ==================== Plot Demand (plot 6) ====================>>

    if ComputePlots == 1
        set(0,'CurrentFigure',1); clf(1); set(gcf,'name','Demand');
        hold on;
        xlabel('Time');
        ylabel('Demand');
        title('Normalized Demand');
        grid minor;
        plot(time,D/DMax);
        plot(time,(D+Export)/DMax);
        legend('Demand','Demand + Exports');
        ylim([0 1]);
        pause(0.001);
    end

    % ==================== Battery ====================>>

    AMax = 1.4e5; % In kWh.
    PAMax = 1e5; % In kW.
    K_A = (PAMax/AMax)/3600; % Normalizing constant.
    A0 = 0.5;
    NAMax = 1;
    NAMin = 0;
    olPn = @(A) min(1,(A-NAMin)/(K_A*dt*10*TMax))*(UseBattery == 1);
    ulPn = @(A) -min(0.35,(NAMax-A)/(K_A*dt*10*TMax))*(UseBattery == 1);

    if UseBattery == 1
        limsA = @(A) [ulPn(A),olPn(A)];
    elseif UseBattery == 0
        limsA = @(A) [];
    else
        disp('Choose UseBattery between 0 or 1.');
        return;
    end

    dA = (NAMax-NAMin)/NA;
    A = NAMin:dA:NAMax;
    if ComputePlots == 1
        set(0,'CurrentFigure',15); clf(15); set(gcf,'name','Battery Control Vs Charge');
        hold on; grid minor;
        xlabel('Capacity');
        ylabel('Control');
        title('Battery Control Limits');
        P = plot(A,olPn(A)); P.LineWidth = 1;
        P = plot(A,ulPn(A)); P.LineWidth = 1;
        legend('Maximum','Minimum');
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
        pause(0.001);
    end

    [vplotA_1,vplot1_A] = meshgrid(A,V1); % For plotting surfaces.
    [vplotA_3,vplot3_A] = meshgrid(A,V3); % For plotting surfaces.
    [vplotA_4,vplot4_A] = meshgrid(A,V4); % For plotting surfaces.

    fA = @(C) TMax*K_A*C;
    KA = 0.3/(MWhTokW); % In USD/kW.

    % ==================== Matrices ====================>>

    u = zeros(length(V1),length(V3),length(V4),length(A));
    CFL = zeros(length(V1),length(V3),length(V4),length(A));
    CFL_v1 = zeros(length(V1),length(V3),length(V4),length(A));
    CFL_v3 = zeros(length(V1),length(V3),length(V4),length(A));
    CFL_v4 = zeros(length(V1),length(V3),length(V4),length(A));
    CFL_A = zeros(length(V1),length(V3),length(V4),length(A));
    HAM = zeros(length(V1),length(V3),length(V4),length(A));
    uV1 = zeros(length(V1),length(V3),length(V4),length(A));
    uV3 = zeros(length(V1),length(V3),length(V4),length(A));
    uV4 = zeros(length(V1),length(V3),length(V4),length(A));
    uA = zeros(length(V1),length(V3),length(V4),length(A));
    ControlsFMC = zeros(length(V1),length(V3),length(V4),length(A),14); % T1,S1,T2,S2,T3,T4,A,F1,F2,F3,F4,V2,V3.
    LV1 = length(V1); % Length of V1.
    LV3 = length(V3); % Length of V3.
    LV4 = length(V4); % Length of V4.
    LA = length(A); % Length of A.
    TL = LV1*LV3*LV4*LA; % Total number of elements in the matrix.
    Controls_T = ControlsFMC;
    UV = {};
    U = {};

    % ==================== More parameters (Water and Fuel) ====================>>

    c1 = -MaxFlowCoeff{4}{1};
    eta1 = MaxFlowCoeff{4}{2};
    d1 = @(V1,tur,spill) c1*(tur*FlowMax1(V1)+ spill*MaxSpill1(V1));
    % Variation in the high due of the released flow (Bonete).

    c2 = -MaxFlowCoeff{4}{3};
    eta2 = MaxFlowCoeff{4}{4};
    d2 = @(V3,tur,spill) c2*(tur*FlowMax2(V3)+spill*MaxSpill2(V3));
    % Variation in the high due of the released flow (Baygorria).

    c3 = -MaxFlowCoeff{4}{5};
    eta3 = MaxFlowCoeff{4}{6};
    d3 = @(V3,tur,spill,t) c3*(tur*FlowMax3(V3,t)+spill*MaxSpill3(V3,t));
    % Variation in the high due of the released flow (Palmar).

    c4 = 0.00147;
    eta4 = 8.82*1.16; % The 1.16 is a correction factor, related with the data of real production.
    d4 = @(tur,spill) c4*(tur*FlowMax4+spill*MaxSpill4);
    % Variation in the high due of the released flow (Salto Grande).

    Khe1 = Matrix{3}(3); % In USD/MWh. Cost per energy of Bonete.
    Khe2 = 0; % In USD/MWh.
    Khe3 = Matrix{3}(2); % In USD/MWh.
    Khe4 = Matrix{3}(1); % In USD/MWh.

    CH1 = Khe1*eta1*(dH1(Vini1)-d1(Vini1,1,0))/(MWhTokW); % Water's value of Bonete.
    CH2 = Khe2*eta2*(dH2(Vini3)-d2(Vini3,1,0))/(MWhTokW); % Water's value of Baygorria.
    CH3 = Khe3*eta3*(dH3(Vini3,1)-d3(Vini3,1,0,1))/(MWhTokW); % Water's value of Palmar.
    CH4 = Khe4*eta4*(dH4(Vini4,1)-d4(1,0))/(MWhTokW); % Water's value of Salto Grande.

    olCT1 = @(uV1,t,V1) FlowMax1(V1)*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCS1 = @(uV1,t,V1) MaxSpill1(V1)*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCT2 = @(uV2,t,V3) FlowMax2(V3)*(CH2 - Lambda32(t+Z32) - uV2/VolMax2 + Lambda21(t));
    olCS2 = @(uV2,t,V3) MaxSpill2(V3)*(CH2 - Lambda32(t+Z32) - uV2/VolMax2 + Lambda21(t));
    olCT3 = @(uV3,t,V3) FlowMax3(V3,t)*(CH3 - uV3/VolMax3);
    olCT4 = @(uV4) FlowMax4*(CH4 - uV4/VolMax4);
    olCA = @(uA) -K_A*uA -KA*1000000*0;
    H0 = @(uV1,uV2,uV3,uV4,t) IT1(t)*uV1/VolMax1 + IT3(t)*uV3/VolMax3 + IT4(t)*uV4/VolMax4 - Lambda21(t)*IT2(t);
    olCV3 = @(uV3,t) MaxTot3*(Lambda32(t) + uV3/VolMax3);
    
    H0_1 = @(uV1,t) IT1(t)*uV1/VolMax1;
    H0_2 = @(uV2,t) - Lambda21(t)*IT2(t);
    H0_3 = @(uV3,t) IT3(t)*uV3/VolMax3;
    H0_4 = @(uV4,t) IT4(t)*uV4/VolMax4;

    k1 = @(t) IT2(t)/MaxV2;
    k2 = @(V3) FlowMax2(V3)/MaxV2;
    k3 = @(V3) MaxSpill2(V3)/MaxV2;
    k = @(t,V3) [k1(t),k2(V3),k3(V3)];

    % Hydraulic power: PH = T^2*K2 + T*K1 + T*S*K3.

    K11 = @(V1) eta1*FlowMax1(V1).*dH1(V1);
    K12 = @(V3) eta2*FlowMax2(V3).*dH2(V3);
    K13 = @(V3,t) eta3*FlowMax3(V3,t).*dH3(V3,t);
    K14 = @(V4,t) eta4*FlowMax4*dH4(V4,t);

    K21 = @(V1) -eta1*FlowMax1(V1).*d1(V1,1,0);
    K22 = @(V3) -eta2*FlowMax2(V3).*d2(V3,1,0);
    K23 = @(V3,t) -eta3*FlowMax3(V3,t).*d3(V3,1,0,t);
    K24 = -eta4*FlowMax4.*d4(1,0);
    
    K31 = @(V1) -eta1*FlowMax1(V1)*d1(V1,0,1);
    K32 = @(V3) -eta2*FlowMax2(V3)*d2(V3,0,1);

    d = @(V1,V3,uV1,uV2,uV3,uV4,uA,t) [olCT1(uV1,t,V1),olCS1(uV1,t,V1),...
        olCT2(uV2,t,V3),olCS2(uV2,t,V3),olCT3(uV3,t,V3),olCT4(uV4),olCA(uA),...
        CostTermicasNorm.*MaxTermicas]';
    dd = diag([0,0,0,0,0,0,KA*1000000*0.000,zeros(1,length(CostTermicasNorm))]);
    rho = KA*1000000*0.07*1*0.1;
    if WhatToDo == 6
        rho = 0;
    end
    b = @(V1,V3,V4,t) [K11(V1),0,K12(V3),0,K13(V3,t),K14(V4,t),PAMax,MaxTermicas]';
    Q = @(V1,V3,t) [K21(V1) K31(V1)/2 0 0 0 0 0 0 0 0 0 0;
        K31(V1)/2 0 0 0 0 0 0 0 0 0 0 0;
        0 0 K22(V3) K32(V3)/2 0 0 0 0 0 0 0 0;
        0 0 K32(V3)/2 0 0 0 0 0 0 0 0 0;
        0 0 0 0 K23(V3,t) 0 0 0 0 0 0 0;
        0 0 0 0 0 K24 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0]; % Is fixed.

    c = @(t) D(t)+Export(t)-WindPower(t)-SolarPower(t)-BioMasPower(t);

    % ==================== Simulation ====================>>
        
    if WhatToDo == 1 || WhatToDo == 3

        tic

        for t = length(time):-1:1

            ControlsFMC_Aux = Controls_T;
            CFL(:) = 0;
            CFL_v1(:) = 0;
            CFL_v3(:) = 0;
            CFL_v4(:) = 0;
            CFL_A(:) = 0;
            
            % ==================== Finite Differences Section ====================>>

            for v1 = 1:length(V1)
                for v3 = 1:length(V3)
                    for v4 = 1:length(V4)
                        for a = 1:length(A)

                            if (v1~=1) && (v3~=1) && (v4~=1) && (a~=1) && (v1~=length(V1)) && (v3~=length(V3)) && (v4~=length(V4)) &&...
                                    (a~=length(A)) && (t<length(time)-1)

                                if (IT1(t+1)-FlowMax1(V1(v1))*ControlsFMC_Aux(v1,v3,v4,a,1)-MaxSpill1(V1(v1))*ControlsFMC_Aux(v1,v3,v4,a,2) >= 0) && Derivatives(1)
                                    uV1(v1,v3,v4,a) = (u(v1+1,v3,v4,a)-u(v1,v3,v4,a))/(dV1);
                                elseif Derivatives(1)
                                    uV1(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1-1,v3,v4,a))/(dV1);
                                end

                                if (IT3(t+1)-FlowMax3(V3(v3),t+1)*ControlsFMC_Aux(v1,v3,v4,a,5)+MaxV3*ControlsFMC_Aux(v1,v3,v4,a,Pc3) >= 0) && Derivatives(2)
                                    uV3(v1,v3,v4,a) = (u(v1,v3+1,v4,a)-u(v1,v3,v4,a))/(dV3);
                                elseif Derivatives(2)
                                    uV3(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1,v3-1,v4,a))/(dV3);
                                end

                                if (IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,a,6) >= 0) && Derivatives(3)
                                    uV4(v1,v3,v4,a) = (u(v1,v3,v4+1,a)-u(v1,v3,v4,a))/(dV4);
                                elseif Derivatives(3)
                                    uV4(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1,v3,v4-1,a))/(dV4);
                                end

                                if (-K_A*ControlsFMC_Aux(v1,v3,v4,a,7) >= 0) && Derivatives(4)
                                    uA(v1,v3,v4,a) = (u(v1,v3,v4,a+1)-u(v1,v3,v4,a))/(dA);
                                elseif Derivatives(4)
                                    uA(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1,v3,v4,a-1))/(dA);
                                end

                            elseif t < length(time) - 1

                                if (v3==1 || v3==length(V3) || v4==1 || v4==length(V4) || a==1 || a==length(A)) && (v1~=1) && (v1~=length(V1)) && Derivatives(1)
                                    if IT1(t+1)-FlowMax1(V1(v1))*ControlsFMC_Aux(v1,v3,v4,a,1)-MaxSpill1(V1(v1))*ControlsFMC_Aux(v1,v3,v4,a,2) >= 0
                                        uV1(v1,v3,v4,a) = (u(v1+1,v3,v4,a)-u(v1,v3,v4,a))/(dV1);
                                    else
                                        uV1(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1-1,v3,v4,a))/(dV1);
                                    end
                                end

                                if (v1==1 || v1==length(V1) || v4==1 || v4==length(V4) || a==1 || a==length(A)) && (v3~=1) && (v3~=length(V3)) && Derivatives(2)
                                    if IT3(t+1)-FlowMax3(V3(v3),t+1)*ControlsFMC_Aux(v1,v3,v4,a,5)+MaxV3*ControlsFMC_Aux(v1,v3,v4,a,Pc3) >= 0
                                        uV3(v1,v3,v4,a) = (u(v1,v3+1,v4,a)-u(v1,v3,v4,a))/(dV3);
                                    else
                                        uV3(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1,v3-1,v4,a))/(dV3);
                                    end
                                end

                                if (v1==1 || v1==length(V1) || v3==1 || v3==length(V3) || a==1 || a==length(A)) && (v4~=1) && (v4~=length(V4)) && Derivatives(3)
                                    if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,a,6) >= 0
                                        uV4(v1,v3,v4,a) = (u(v1,v3,v4+1,a)-u(v1,v3,v4,a))/(dV4);
                                    else
                                        uV4(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1,v3,v4-1,a))/(dV4);
                                    end
                                end

                                if (v1==1 || v1==length(V1) || v3==1 || v3==length(V3) || v4==1 || v4==length(V4)) && (a~=1) && (a~=length(A)) && Derivatives(4)
                                    if -K_A*ControlsFMC_Aux(v1,v3,v4,a,7) >= 0
                                        uA(v1,v3,v4,a) = (u(v1,v3,v4,a+1)-u(v1,v3,v4,a))/(dA);
                                    else
                                        uA(v1,v3,v4,a) = (u(v1,v3,v4,a)-u(v1,v3,v4,a-1))/(dA);
                                    end
                                end

                            end

                            CFL(v1,v3,v4,a) = 1 / (abs(f1(t,V1(v1),ControlsFMC_Aux(v1,v3,v4,a,1),ControlsFMC_Aux(v1,v3,v4,a,2)))/dV1 + ...
                                abs(f3(t,V3(v3),ControlsFMC_Aux(v1,v3,v4,a,5),0,ControlsFMC_Aux(v1,v3,v4,a,Pc3)))/dV3 + ...
                                abs(f4(t,V4(v4),ControlsFMC_Aux(v1,v3,v4,a,6),0))/dV4 + ...
                                abs(fA(ControlsFMC_Aux(v1,v3,v4,a,7)))/dA);
                            CFL_v1(v1,v3,v4,a) = abs(f1(t,V1(v1),ControlsFMC_Aux(v1,v3,v4,a,1),ControlsFMC_Aux(v1,v3,v4,a,2)))/dV1;
                            CFL_v3(v1,v3,v4,a) = abs(f3(t,V3(v3),ControlsFMC_Aux(v1,v3,v4,a,5),0,ControlsFMC_Aux(v1,v3,v4,a,Pc3)))/dV3;
                            CFL_v4(v1,v3,v4,a) = abs(f4(t,V4(v4),ControlsFMC_Aux(v1,v3,v4,a,6),0))/dV4;
                            CFL_A(v1,v3,v4,a) = abs(fA(ControlsFMC_Aux(v1,v3,v4,a,7)))/dA;

                        end
                    end
                end
            end

            uV1(1,:,:,:) = uV1(2,:,:,:);
            uV1(end,:,:,:) = uV1(end-1,:,:,:);
            uV3(:,1,:,:) = uV3(:,2,:,:);
            uV3(:,end,:,:) = uV3(:,end-1,:,:);
            uV4(:,:,1,:) = uV4(:,:,2,:);
            uV4(:,:,end,:) = uV4(:,:,end-1,:);
            uA(:,:,:,1) = uA(:,:,:,2);
            uA(:,:,:,end) = uA(:,:,:,end-1);
            UV{t} = {uV1,uV3,uV4,uA};
            
            % ==================== Finite Differences Section ====================>>

            fprintf('Delta T / CFL = %.8f.\n',dt/min(min(min(min(CFL)))));
            
            CFL_Plot(t,:) = [dt/min(min(min(min(CFL)))),max(max(max(max(CFL_v1)))),...
                max(max(max(max(CFL_v3)))),max(max(max(max(CFL_v4)))),max(max(max(max(CFL_A))))];
            
            % ==================== Controls Section ====================>>

            Size_Mat = size(u);
            uV1_PF = uV1(:);
            uV3_PF = uV3(:);
            uV4_PF = uV4(:);
            Controls_Aux = zeros(TL,14);
            uA_PF = uA(:);
            Warm_Start = reshape(Controls_T,LV1,LV3,LV4,LA,14);
            Warm_Start = Warm_Start(:,:,:,:,1:12);

            if WarmStart == 0
                Warm_Start(Warm_Start ~= 0) = 0;
            elseif WarmStart ~= 1
                disp('Choose WarmStart 0 or 1.');
                return;
            end

            parfor index = 1:TL

                v3_aux = uV3_PF(index);
                
                if (UseLamb32 == 0) || (olCV3(uV3(index),t) >= 0 && t > Z32)
                    Aux_V3 = 0;
                elseif t > Z32
                    if v3_aux > 0
                        Aux_V3 = 1;
                    else
                        Aux_V3 = 0;
                    end
                else
                    Aux_V3 = Phi_V3(t);
                end

                if UseLamb21 == 1
                    if t > Z21
                        Aux_V2 = 0;
                        PhiV2Fixed = 0;
                    else
                        Aux_V2 = Phi_V2(t);
                        PhiV2Fixed = 1;
                    end
                else
                    Aux_V2 = 0;
                    PhiV2Fixed = 1;
                end

                [ind_V1,ind_V3,ind_V4,ind_A] = ind2sub(Size_Mat,index);
                WS = squeeze(Warm_Start(ind_V1,ind_V3,ind_V4,ind_A,:));
                Aux_uA = uA_PF(index);
                Aux_Lims = limsA(A(ind_A));

                Controls_Aux(index,:) = [Quad_FMC_Lambda_7(Q(V1(ind_V1),V3(ind_V3),t),b(V1(ind_V1),V3(ind_V3),...
                        V4(ind_V4),t),c(t),d(V1(ind_V1),V3(ind_V3),uV1_PF(index),0,v3_aux,uV4_PF(index),...
                        Aux_uA,t),dd,k(t,V3(ind_V3)),options,PhiV2Fixed,Aux_V2,WS,Aux_Lims,InfTermicas);Aux_V3];

            end

            Controls_T = reshape(Controls_Aux,LV1,LV3,LV4,LA,14);

            % ==================== Time Step ====================>>

            if t ~= length(time)
                for v1 = 1:length(V1)
                    for v3 = 1:length(V3)
                        for v4 = 1:length(V4)
                            for a = 1:length(A)

                                Controls(1:14) = ControlsFMC_Aux(v1,v3,v4,a,:);

                                HAM(v1,v3,v4,a) = TMax*(Controls*[d(V1(v1),V3(v3),uV1(v1,v3,v4,a),...
                                    0,uV3(v1,v3,v4,a),uV4(v1,v3,v4,a),uA(v1,v3,v4,a),t);0;...
                                    olCV3(uV3(v1,v3,v4,a),t)]+H0(uV1(v1,v3,v4,a),0,uV3(v1,v3,v4,a),uV4(v1,v3,v4,a),t));

                                u(v1,v3,v4,a) = u(v1,v3,v4,a) + dt*HAM(v1,v3,v4,a);

                            end
                        end
                    end
                end
            end
            
            U{t} = u;

            % ==================== Plotting over time ====================>>

            if ComputePlots == 1

                set(0,'CurrentFigure',1);
                vline(time(t),'r',' ');
                title(['Normalized Demand (time=',num2str(t),')']);      
                
                set(0,'CurrentFigure',2);
                clf(2); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(u(index,:,:,1)/NormCost)');
                    title(['Optimal Cost (A = 0) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',3);
                clf(3); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(u(index,:,:,end)/NormCost)');
                    title(['Optimal Cost (A = 1) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',4);
                clf(4); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(u(index,:,:,1+floor(length(A)/2))/NormCost)');
                    title(['Optimal Cost (A = 0.5) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',5);
                clf(5); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(HAM(index,:,:,1)/NormCost)');
                    title(['Hamiltonian (A = 0) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',6);
                clf(6); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(HAM(index,:,:,end)/NormCost)');
                    title(['Hamiltonian (A = 1) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',7);
                clf(7); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(HAM(index,:,:,1+floor(length(A)/2))/NormCost)');
                    title(['Hamiltonian (A = 0.5) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',8);
                clf(8); hold on;
                for i = 1:floor(length(V3)/2):length(V3)
                    surf(vplot4_A,vplotA_4,squeeze(u(end,i,:,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['Optimal Cost ($V1 = \underline{V}_1$) at time = ',num2str(t)],'Interpreter','latex');
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',9);
                clf(9); hold on;
                for index = 1:floor(length(V4)/2):length(V4)
                    surf(vplot3_A,vplotA_3,squeeze(uV3(index,:,index,:)/NormCost));
                end
                xlabel('Vol. Palmar'); ylabel('Charge Battery');
                title(['u_{V3} at time = ',num2str(t)]);
                xlim([Vmin3 Vmax3]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',10);
                clf(10); hold on;
                for index = 1:length(V4)
                    surf(vplot3_A,vplotA_3,squeeze(uV1(index,:,index,:)/NormCost));
                end
                xlabel('Vol. Palmar'); ylabel('Charge Battery');
                title(['u_{V1} at time = ',num2str(t)]);
                xlim([Vmin3 Vmax3]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',14);
                clf(14); hold on;
                for i = 1:floor(length(V4)/2):length(V4)
                    surf(vplot4_A,vplotA_4,squeeze(u(i,i,:,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['Optimal Cost at time = ',num2str(t)]);
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',12);
                clf(12); hold on;
                for index = 1:length(V4)
                    surf(vplot4_A,vplotA_4,squeeze(uA(index,:,index,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['u_{A} at time = ',num2str(t)]);
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',13);
                clf(13); hold on;
                Vols = linspace(min([Vmin1,Vmin3,Vmin4]),max([Vmax1,Vmax3,Vmax4]),1000);
                plot(Vols,K11(Vols)+K21(Vols),Vols,K12(Vols)+K22(Vols),...
                    Vols,K13(Vols,t)+K23(Vols,t),Vols,K14(Vols,t)+K24);
                P = plot(Vols,K12(Vols)+K22(Vols)+K13(Vols,t)+K23(Vols,t));
                P.LineWidth = 1;
                xlabel('Volume'); ylabel('Power (kW)');
                legend('Power Bonete','Power Baygorria','Power Palmar','Power Salto Grande',...
                    'Power Baygorria + Palmar');
                title(['Power at time = ',num2str(t)]);
                xlim([min(Vols) max(Vols)]);

            end

            fprintf('Discretizations = %d. Completed: %.2f.\n',NV4,100*(1-t/length(time)));
            pause(0.1);

        end

        IP1 = min(abs(V1-Vini1)) == abs(V1-Vini1);
        IP3 = min(abs(V3-Vini3)) == abs(V3-Vini3);
        IP4 = min(abs(V4-Vini4)) == abs(V4-Vini4);
        IPA = min(abs(A-A0)) == abs(A-A0);
        Final_Cost = u(IP1,IP3,IP4,IPA)/NormCost;

        ToSave = {Final_Cost,UV,CFL_Plot,U};
        save([pwd '/','Simulation_',name,'.mat'],'ToSave');
        toc
        
    end
    
    %% Optimal Path Section:
    
    if WhatToDo == 2 || WhatToDo == 3 || WhatToDo == 4 || WhatToDo == 5 || WhatToDo == 6 || WhatToDo == 7
        
        tic
        
        if WhatToDo == 6
            Fixed = 0;
        else
            Fixed = 1;
        end
            
        load([pwd '/','Simulation_',name,'.mat']);
        if WhatToDo == 4
            load([pwd '/','Historical/Energy_Controls_',num2str(Day),'_2019.mat']);
        end
        if WhatToDo == 7
            load([pwd '/','Simulation_Controls_',name,'_Day_',num2str(Day),'.mat']);
        end
        
        Final_Cost = ToSave{1};
        UV = ToSave{2};
        CFL_Plot = ToSave{3};
        
        if WhatToDo ~= 5
        
            IP1 = min(abs(V1-Vini1)) == abs(V1-Vini1);
            IP3 = min(abs(V3-Vini3)) == abs(V3-Vini3);
            IP4 = min(abs(V4-Vini4)) == abs(V4-Vini4);
            IPA = min(abs(A-A0)) == abs(A-A0);

            Optimal_Cost = zeros(1,length(time));
            Phi_V2_Aux = zeros(1,length(time)); % To save the past virtual control V2.
            Phi_V3_Aux = zeros(1,length(time)); % To save the past virtual control V3.
            Power_H1_OP = zeros(1,length(time));
            Power_H2_OP = zeros(1,length(time));
            Power_H3_OP = zeros(1,length(time));
            Power_H4_OP = zeros(1,length(time));
            Vol_OP_1 = zeros(1,length(time));
            Vol_OP_3 = zeros(1,length(time));
            Vol_OP_4 = zeros(1,length(time));
            Vol_OP_Real_1 = zeros(1,length(time));
            Vol_OP_Real_3 = zeros(1,length(time));
            Vol_OP_Real_4 = zeros(1,length(time));
            Vol_OP_1(1) = V1(IP1);
            Vol_OP_3(1) = V3(IP3);
            Vol_OP_4(1) = V4(IP4);
            Vol_OP_Real_1(1) = V1(IP1);
            Vol_OP_Real_3(1) = V3(IP3);
            Vol_OP_Real_4(1) = V4(IP4);
            Cost = zeros(10,length(time));
            uV1plot = zeros(1,length(time));
            uV3plot = zeros(1,length(time));
            uV4plot = zeros(1,length(time));
            CostAcum1 = 0;
            CostAcum2 = 0;
            CostAcum3 = 0;
            CostAcum4 = 0;
            CostAcum5 = 0;
            CostAcum6 = 0;
            CostAcum7 = 0;
            CostAcum8 = 0;
            CostAcum9 = 0;
            CostAcum1R = 0;
            CostAcum2R = 0;
            CostAcum3R = 0;
            CostAcum4R = 0;
            CostAcum5R = 0;
            CostAcum6R = 0;
            CostAcum7R = 0;
            CostAcum8R = 0;
            CostAcum9R = 0;
            NC_Cost_Acum = 0;
            if WhatToDo ~= 7
                Controls_OP = zeros(length(time),14);
            end
            CostControls = zeros(12,length(time));
            CostControlsInstant = zeros(12,length(time));
            Battery_Real = zeros(1,length(time));
            Power_Battery_OP = zeros(1,length(time));
            Prev_SG = zeros(1,length(time)+1);
            Prev_Batt = zeros(1,length(time)+1);
            Prev_F4 = zeros(1,length(time)+1);
            Prec_Cont = zeros(length(time)+1,12);
            Power_F1_OP = zeros(1,length(time));
            Power_F2_OP = zeros(1,length(time));
            Power_F3_OP = zeros(1,length(time));
            Power_F4_OP = zeros(1,length(time));
            Power_F5_OP = zeros(1,length(time));
            Battery_Cost = 0;
            Extra_Cost = 0;
            lambda = {};

            uAplot = zeros(1,length(time));
            Battery_OP = A(IPA);
            Battery_Real(1) = A(IPA);
            x0 = zeros(12,1);

            Power_H1 = @(x1,x2,V1) eta1*x1*FlowMax1(V1)*(dH1(V1)-d1(V1,x1,x2)); % In kW.
            Power_H2 = @(x3,x4,V3) eta2*x3*FlowMax2(V3)*(dH2(V3)-d2(V3,x3,x4)); % In kW.
            Power_H3 = @(x5,V3,t) eta3*x5*FlowMax3(V3,t)*(dH3(V3,t)-d3(V3,x5,0,t)); % In kW.
            Power_H4 = @(x6,V4,t) eta4*x6*FlowMax4*(dH4(V4,t)-d4(x6,0)); % In kW.
            Power_Battery = @(x7) x7*PAMax; % In kW.

            for t = 1:length(time)

                fprintf('Completed: %.2f.\n',100*(t/length(time)));

                if t == 1
                    uV1 = UV{t}{1}(IP1,IP3,IP4,IPA);
                    uV3 = UV{t}{2}(IP1,IP3,IP4,IPA);
                    uV4 = UV{t}{3}(IP1,IP3,IP4,IPA);
                    uA = UV{t}{4}(IP1,IP3,IP4,IPA);
                else
                    realxV1 = (Vol_OP_Real_1(t)-V1(IP1))/(V1(IP1+1)-V1(IP1));
                    realxV3 = (Vol_OP_Real_3(t)-V3(IP3))/(V3(IP3+1)-V3(IP3));
                    realxV4 = (Vol_OP_Real_4(t)-V4(IP4))/(V4(IP4+1)-V4(IP4));
                    realxA = (Battery_Real(t)-A(IPA))/(A(IPA+1)-A(IPA));
                    real = [realxV1,realxV3,realxV4,realxA];

                    valsV1 = [UV{t}{1}(IP1,IP3,IP4,IPA);UV{t}{1}(IP1+1,IP3,IP4,IPA);UV{t}{1}(IP1,IP3+1,IP4,IPA);...
                        UV{t}{1}(IP1+1,IP3+1,IP4,IPA);UV{t}{1}(IP1,IP3,IP4+1,IPA);UV{t}{1}(IP1+1,IP3,IP4+1,IPA);...
                        UV{t}{1}(IP1,IP3+1,IP4+1,IPA);UV{t}{1}(IP1+1,IP3+1,IP4+1,IPA);UV{t}{1}(IP1,IP3,IP4,IPA+1);...
                        UV{t}{1}(IP1+1,IP3,IP4,IPA+1);UV{t}{1}(IP1,IP3+1,IP4,IPA+1);UV{t}{1}(IP1+1,IP3+1,IP4,IPA+1);...
                        UV{t}{1}(IP1,IP3,IP4+1,IPA+1);UV{t}{1}(IP1+1,IP3,IP4+1,IPA+1);UV{t}{1}(IP1,IP3+1,IP4+1,IPA+1);...
                        UV{t}{1}(IP1+1,IP3+1,IP4+1,IPA+1)];
                    uV1 = Int4D_2(real,valsV1);

                    valsV3 = [UV{t}{2}(IP1,IP3,IP4,IPA);UV{t}{2}(IP1+1,IP3,IP4,IPA);UV{t}{2}(IP1,IP3+1,IP4,IPA);...
                        UV{t}{2}(IP1+1,IP3+1,IP4,IPA);UV{t}{2}(IP1,IP3,IP4+1,IPA);UV{t}{2}(IP1+1,IP3,IP4+1,IPA);...
                        UV{t}{2}(IP1,IP3+1,IP4+1,IPA);UV{t}{2}(IP1+1,IP3+1,IP4+1,IPA);UV{t}{2}(IP1,IP3,IP4,IPA+1);...
                        UV{t}{2}(IP1+1,IP3,IP4,IPA+1);UV{t}{2}(IP1,IP3+1,IP4,IPA+1);UV{t}{2}(IP1+1,IP3+1,IP4,IPA+1);...
                        UV{t}{2}(IP1,IP3,IP4+1,IPA+1);UV{t}{2}(IP1+1,IP3,IP4+1,IPA+1);UV{t}{2}(IP1,IP3+1,IP4+1,IPA+1);...
                        UV{t}{2}(IP1+1,IP3+1,IP4+1,IPA+1)];
                    uV3 = Int4D_2(real,valsV3);

                    valsV4 = [UV{t}{3}(IP1,IP3,IP4,IPA);UV{t}{3}(IP1+1,IP3,IP4,IPA);UV{t}{3}(IP1,IP3+1,IP4,IPA);...
                        UV{t}{3}(IP1+1,IP3+1,IP4,IPA);UV{t}{3}(IP1,IP3,IP4+1,IPA);UV{t}{3}(IP1+1,IP3,IP4+1,IPA);...
                        UV{t}{3}(IP1,IP3+1,IP4+1,IPA);UV{t}{3}(IP1+1,IP3+1,IP4+1,IPA);UV{t}{3}(IP1,IP3,IP4,IPA+1);...
                        UV{t}{3}(IP1+1,IP3,IP4,IPA+1);UV{t}{3}(IP1,IP3+1,IP4,IPA+1);UV{t}{3}(IP1+1,IP3+1,IP4,IPA+1);...
                        UV{t}{3}(IP1,IP3,IP4+1,IPA+1);UV{t}{3}(IP1+1,IP3,IP4+1,IPA+1);UV{t}{3}(IP1,IP3+1,IP4+1,IPA+1);...
                        UV{t}{3}(IP1+1,IP3+1,IP4+1,IPA+1)];
                    uV4 = Int4D_2(real,valsV4);

                    valsA = [UV{t}{4}(IP1,IP3,IP4,IPA);UV{t}{4}(IP1+1,IP3,IP4,IPA);UV{t}{4}(IP1,IP3+1,IP4,IPA);...
                        UV{t}{4}(IP1+1,IP3+1,IP4,IPA);UV{t}{4}(IP1,IP3,IP4+1,IPA);UV{t}{4}(IP1+1,IP3,IP4+1,IPA);...
                        UV{t}{4}(IP1,IP3+1,IP4+1,IPA);UV{t}{4}(IP1+1,IP3+1,IP4+1,IPA);UV{t}{4}(IP1,IP3,IP4,IPA+1);...
                        UV{t}{4}(IP1+1,IP3,IP4,IPA+1);UV{t}{4}(IP1,IP3+1,IP4,IPA+1);UV{t}{4}(IP1+1,IP3+1,IP4,IPA+1);...
                        UV{t}{4}(IP1,IP3,IP4+1,IPA+1);UV{t}{4}(IP1+1,IP3,IP4+1,IPA+1);UV{t}{4}(IP1,IP3+1,IP4+1,IPA+1);...
                        UV{t}{4}(IP1+1,IP3+1,IP4+1,IPA+1)];
                    uA = Int4D_2(real,valsA);

                end

                ContLims = limsA(Battery_Real(t));
                uV2 = 0;
                uV1plot(t) = uV1;
                uV3plot(t) = uV3;
                uV4plot(t) = uV4;
                uAplot(t) = uA;

                % ==================== Controls Section ====================>>

                if t > Z21
                    if UseLamb21 == 1
                        Real2 = Phi_V2_Aux(t-Z21);
                    else
                        Real2 = 0;
                    end
                else
                    if UseLamb21 == 1
                        Real2 = Phi_V2(t);
                    else
                        Real2 = 0;
                    end
                end

                if t == 1
                    dd_OP = dd;
                    d_OP = @(V1,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V3,uV1,uV2,uV3,uV4,uA,t);
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t);
                else
                    dd_OP = dd + diag(rho/dt*ones(1,12)) - diag(rho/dt*ones(1,12).*[0,1,0,0,0,0,0,0,0,0,0,0])*9900/10000;
                    d_OP = @(V1,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V3,uV1,uV2,uV3,uV4,uA,t) + ...
                    -2*rho/dt*Conts' + 2*rho/dt*Conts(2)*[0,1,0,0,0,0,0,0,0,0,0,0]';
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t) + rho*sum(Conts.^2)/dt - rho*sum(Conts(2)^2)/dt*[0,1,0,0,0,0,0,0,0,0,0,0]';
                    Added_Cost = @(New_Conts,Conts) sum(New_Conts.^2)*rho/dt - 2*rho/dt*sum(New_Conts.*Conts) + rho/dt*sum(Conts.^2); % <<< I have to correct that, the spillage of Bonete now is L2.
                end

                if WhatToDo ~= 4 && WhatToDo ~= 7

                    [Controls_OP(t,1:13)] = Quad_FMC_Lambda_7(Q(Vol_OP_Real_1(t),Vol_OP_Real_3(t),t),...
                        b(Vol_OP_Real_1(t),Vol_OP_Real_3(t),Vol_OP_Real_4(t),t),c(t),...
                        d_OP(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t,Prec_Cont(t,:)),dd_OP,k(t,Vol_OP_Real_3(t)),...
                        options,Fixed,Real2,x0,ContLims,InfTermicas);

                    dd_L = dd;
                    d_L = @(V1,V3,uV1,uV2,uV3,uV4,uA,t,ContSG,ContA,ContF4) d(V1,V3,uV1,uV2,uV3,uV4,uA,t);
                    [~,lambda{t}] = Quad_FMC_Lambda_6(Q(Vol_OP_Real_1(t),Vol_OP_Real_3(t),t),...
                        b(Vol_OP_Real_1(t),Vol_OP_Real_3(t),Vol_OP_Real_4(t),t),c(t),...
                        d_L(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t,Prec_Cont(t,:)),dd_L,k(t,Vol_OP_Real_3(t)),...
                        options,Fixed,Real2,Controls_OP(t,1:12)',ContLims,InfTermicas);

                elseif WhatToDo == 4

                    Controls_OP(t,1:13) = [Save{2}{1}(t),0,Save{2}{2}(t),0,Save{2}{3}(t),Save{2}{4}(t),0,Save{2}{5}(t),...
                        Save{2}{6}(t),Save{2}{7}(t),Save{2}{8}(t),0,0];
                    Controls_OP(t,1:12) = Controls_OP(t,1:12)*1000./[Power_H1(1,0,Vol_OP_Real_1(t)),1,...
                    Power_H2(1,0,Vol_OP_Real_3(t)),1,...
                    Power_H3(1,Vol_OP_Real_3(t),t),Power_H4(1,Vol_OP_Real_4(t),t),PAMax,MaxTermicas];

                elseif WhatToDo == 7
                    
                end

                Phi_V2_Aux(t) = (Controls_OP(t,1)*FlowMax1(Vol_OP_Real_1(t)) + Controls_OP(t,2)*MaxSpill1(Vol_OP_Real_1(t))) / MaxV2;
                Phi_V3_Aux(t) = (Controls_OP(t,3)*FlowMax2(Vol_OP_Real_3(t)) + Controls_OP(t,4)*MaxSpill2(Vol_OP_Real_3(t))) / MaxV3;
                x0(1:12) = Controls_OP(t,1:12);
                Prec_Cont(t+1,1:12) = Controls_OP(t,1:12);

                if UseLamb32 == 1
                    if t > Z32
                        Controls_OP(t,14) = Phi_V3_Aux(t-Z32);
                        if Fixed == 0
                            if olCV3(uV3,t) > 0
                                Controls_OP(t,14) = 0;
                            else
                                Controls_OP(t,14) = 1;
                            end
                        end
                    else
                        Controls_OP(t,14) = Phi_V3(t);
                    end
                else
                    Controls_OP(t,14) = 0;
                end

                WaterBonete(t) =  dt*TMax*(Controls_OP(t,1)*FlowMax1(Vol_OP_Real_1(t)) + Controls_OP(t,2)*MaxSpill1(Vol_OP_Real_1(t))); 
                WaterVirtualBay(t) = dt*TMax*Controls_OP(t,13)*MaxV2;

                Aux_Cost = dt*TMax*Controls_OP(t,1:12).*d(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t)';
                NC_Cost(t) = dt*TMax*H0(uV1,uV2,uV3,uV4,t);

                Cost(1,t) = Aux_Cost(1)+Aux_Cost(2);
                Cost(2,t) = Aux_Cost(3)+Aux_Cost(4);
                Cost(3,t) = Aux_Cost(5);
                Cost(4,t) = Aux_Cost(6);
                Cost(5,t) = Aux_Cost(8);
                Cost(6,t) = Aux_Cost(9);
                Cost(7,t) = Aux_Cost(10);
                Cost(8,t) = Aux_Cost(11);
                Cost(9,t) = Aux_Cost(7);
                Cost(10,t) = Aux_Cost(12);
                
                CostR(1,t) = Cost(1,t) + TMax*dt*H0_1(uV1,t);
                CostR(2,t) = Cost(2,t) + TMax*dt*H0_2(uV2,t);
                CostR(3,t) = Cost(3,t) + TMax*dt*H0_3(uV3,t);
                CostR(4,t) = Cost(4,t) + TMax*dt*H0_4(uV4,t);
                
                if t > 1
                    CostAcum1(t) = CostAcum1(t-1) + Cost(1,t);
                    CostAcum2(t) = CostAcum2(t-1) + Cost(2,t);
                    CostAcum3(t) = CostAcum3(t-1) + Cost(3,t);
                    CostAcum4(t) = CostAcum4(t-1) + Cost(4,t);
                    CostAcum5(t) = CostAcum5(t-1) + Cost(5,t);
                    CostAcum6(t) = CostAcum6(t-1) + Cost(6,t);
                    CostAcum7(t) = CostAcum7(t-1) + Cost(7,t);
                    CostAcum8(t) = CostAcum8(t-1) + Cost(8,t);
                    CostAcum9(t) = CostAcum9(t-1) + Cost(10,t);
                    NC_Cost_Acum(t) = NC_Cost_Acum(t-1) + NC_Cost(t);
                    CostAcum1R(t) = CostAcum1R(t-1) + CostR(1,t);
                    CostAcum2R(t) = CostAcum2R(t-1) + CostR(2,t);
                    CostAcum3R(t) = CostAcum3R(t-1) + CostR(3,t);
                    CostAcum4R(t) = CostAcum4R(t-1) + CostR(4,t);
                    CostAcum5R(t) = CostAcum5R(t-1) + Cost(5,t);
                    CostAcum6R(t) = CostAcum6R(t-1) + Cost(6,t);
                    CostAcum7R(t) = CostAcum7R(t-1) + Cost(7,t);
                    CostAcum8R(t) = CostAcum8R(t-1) + Cost(8,t);
                    CostAcum9R(t) = CostAcum9R(t-1) + Cost(10,t);
                end
                
                CostControls(:,t) = d(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t)./...
                   [Power_H1(1,Controls_OP(t,2),Vol_OP_Real_1(t)),MaxSpill1(Vol_OP_Real_1(t))*dt*TMax/m3Tohm3,...
                   Power_H2(1,Controls_OP(t,4),Vol_OP_Real_3(t)),MaxSpill2(Vol_OP_Real_3(t))*dt*TMax/m3Tohm3,...
                   Power_H3(1,Vol_OP_Real_3(t),t),Power_H4(1,Vol_OP_Real_4(t),t),PAMax,MaxTermicas]';
               % This cost is in USD/kW for all but the spillages, and is the cost when the devices are working
               % at their maximum power. For the spillages is USD/Control.

                CostControlsInstant(:,t) = d(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t).*Controls_OP(t,1:12)'./...
                   ([Power_H1(Controls_OP(t,1),Controls_OP(t,2),Vol_OP_Real_1(t)),MaxSpill1(Vol_OP_Real_1(t))*dt*TMax/m3Tohm3,...
                   Power_H2(Controls_OP(t,3),Controls_OP(t,4),Vol_OP_Real_3(t)),MaxSpill2(Vol_OP_Real_3(t))*dt*TMax/m3Tohm3,...
                   Power_H3(Controls_OP(t,5),Vol_OP_Real_3(t),t),Power_H4(Controls_OP(t,6),Vol_OP_Real_4(t),t),PAMax,MaxTermicas]'+...
                   ones(12,1));
               % This cost is in USD/kW for the optimal path. We add at the end
               % ones(11,1) to avoid 0/0 when come control is = 0.

                % ==================== Dynamics and Powers Section ====================>>

                if t ~= length(time)

                    Vol_OP_Real_1(t+1) = Vol_OP_Real_1(t) + (dt*TMax/VolMax1)*(IT1(t)-...
                        Controls_OP(t,1)*FlowMax1(Vol_OP_1(t))-Controls_OP(t,2)*MaxSpill1(Vol_OP_1(t)));
                    Vol_OP_Real_3(t+1) = Vol_OP_Real_3(t) + (dt*TMax/VolMax3)*(IT3(t)-...
                        Controls_OP(t,5)*FlowMax3(Vol_OP_3(t),t)+Controls_OP(t,14)*MaxV3);
                    Vol_OP_Real_4(t+1) = Vol_OP_Real_4(t) + (dt*TMax/VolMax4)*(IT4(t)-...
                        Controls_OP(t,6)*FlowMax4);
                    Battery_Real(t+1) = Battery_Real(t) - dt*TMax*K_A*Controls_OP(t,7);

                    AuxV1 = V1(min(abs(V1-Vol_OP_Real_1(t+1))) == abs(V1-Vol_OP_Real_1(t+1)));
                    Vol_OP_1(t+1) = AuxV1(1);
                    AuxV3 = V3(min(abs(V3-Vol_OP_Real_3(t+1))) == abs(V3-Vol_OP_Real_3(t+1)));
                    Vol_OP_3(t+1) = AuxV3(1);
                    AuxV4 = V4(min(abs(V4-Vol_OP_Real_4(t+1))) == abs(V4-Vol_OP_Real_4(t+1)));
                    Vol_OP_4(t+1) = AuxV4(1);
                    AuxA = A(min(abs(A-Battery_Real(t+1))) == abs(A-Battery_Real(t+1)));
                    Battery_OP(t+1) = AuxA(1);

                    IP1 = find(V1 == Vol_OP_1(t+1));
                    IP3 = find(V3 == Vol_OP_3(t+1));
                    IP4 = find(V4 == Vol_OP_4(t+1));
                    IPA = find(A == Battery_OP(t+1));

                    if Vol_OP_1(t+1) > Vol_OP_Real_1(t+1)
                        IP1 = IP1 - 1;
                    end
                    if Vol_OP_3(t+1) > Vol_OP_Real_3(t+1)
                        IP3 = IP3 - 1;
                    end
                    if Vol_OP_4(t+1) > Vol_OP_Real_4(t+1)
                        IP4 = IP4 - 1;
                    end
                    if Battery_OP(t+1) > Battery_Real(t+1)
                        IPA = IPA - 1;
                    end

                    if Battery_Real(t+1) >= NAMax
                        Battery_Real(t+1) = NAMax;
                        IPA = length(A)-1;
                    elseif Battery_Real(t+1) <= NAMin
                        Battery_Real(t+1) = NAMin;
                        IPA = 1;
                    end

                    if Vol_OP_Real_1(t+1) >= Vmax1
                        Vol_OP_Real_1(t+1) = Vmax1;
                        IP1 = length(V1)-1;
                    elseif Vol_OP_Real_1(t+1) <= Vmin1
                        Vol_OP_Real_1(t+1) = Vmin1;
                        IP1 = 1;
                    end

                    if Vol_OP_Real_3(t+1) >= Vmax3
                        Vol_OP_Real_3(t+1) = Vmax3;
                        IP3 = length(V3)-1;
                    elseif Vol_OP_Real_3(t+1) <= Vmin3
                        Vol_OP_Real_3(t+1) = Vmin3;
                        IP3 = 1;
                    end

                    if Vol_OP_Real_4(t+1) >= Vmax4
                        Vol_OP_Real_4(t+1) = Vmax4;
                        IP4 = length(V4)-1;
                    elseif Vol_OP_Real_4(t+1) <= Vmin4
                        Vol_OP_Real_4(t+1) = Vmin4;
                        IP4 = 1;
                    end

                    dd_OP = dd;
                    d_OP = @(V1,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V3,uV1,uV2,uV3,uV4,uA,t);
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t);

                    Optimal_Cost(t+1) = Optimal_Cost(t) + dt*TMax*(Controls_OP(t,:)*...
                        [d_OP(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t,Prec_Cont(t,:));0;olCV3(uV3,t)] +...
                        Controls_OP(t,1:12)*dd_OP*Controls_OP(t,1:12)' +...
                        H0_OP(uV1,uV2,uV3,uV4,t,Prec_Cont(t,:)));

                    Battery_Cost = Battery_Cost + dt*TMax*Controls_OP(t,7)*olCA(uA);

                end

                Power_H1_OP(t) = Power_H1(Controls_OP(t,1),Controls_OP(t,2),Vol_OP_Real_1(t));
                Power_H2_OP(t) = Power_H2(Controls_OP(t,3),Controls_OP(t,4),Vol_OP_Real_3(t));
                Power_H3_OP(t) = Power_H3(Controls_OP(t,5),Vol_OP_Real_3(t),t);
                Power_H4_OP(t) = Power_H4(Controls_OP(t,6),Vol_OP_Real_4(t),t);
                Power_Battery_OP(t) = Power_Battery(Controls_OP(t,7));
                Power_F1_OP(t) = MaxTermicas(1)*(Controls_OP(t,8));
                Power_F2_OP(t) = MaxTermicas(2)*(Controls_OP(t,9));
                Power_F3_OP(t) = MaxTermicas(3)*(Controls_OP(t,10));
                Power_F4_OP(t) = MaxTermicas(4)*(Controls_OP(t,11));
                Power_F5_OP(t) = MaxTermicas(5)*(Controls_OP(t,12));

            end

            % ==================== 'Others' Section ====================>>

            if ComputePlots == 1

                if WhatToDo ~= 4 && WhatToDo ~= 6 && WhatToDo ~= 7
                    for i = 1:length(time)
                        Spot(i) = lambda{i}.ineqnonlin*MWhTokW;
                    end
                end
                
                Format_Long = 1;
                
                if Format_Long == 0
                    set(0,'CurrentFigure',101); clf(101); set(gcf,'Position',[10 10 900 420],'name','Power Balance');
                else
                	set(0,'CurrentFigure',101); clf(101); set(gcf,'Position',[10 10 750*2 300*2],'name','Power Balance');
                end
                    grid minor; hold on;
                if UseBattery == 1
                    hh = area(time,[BioMasPower;SolarPower;WindPower;Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                    Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP;Power_Battery_OP.*(Power_Battery_OP > 0)]'/1000);
                    hh(13).FaceColor = [1 1 0]; % Battery.
                else
                    hh = area(time,[BioMasPower;SolarPower;WindPower;Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                    Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP]'/1000);
                end
                hh(1).FaceColor = [1 .5 0]; % Biomasa.
                hh(2).FaceColor = [1 1 .5]; % Solar.
                hh(3).FaceColor = [0 1 0]; % Wind.
                hh(4).FaceColor = [.5 .5 1]; % Bonete.
                hh(5).FaceColor = [.5 .8 1]; % Baygorria.
                hh(6).FaceColor = [.2 .2 1]; % Palmar.
                hh(7).FaceColor = [.8 .8 1]; % SG.
                hh(8).FaceColor = [1 .6 .6]; % T1.
                hh(9).FaceColor = [1 .3 .5]; % T2.
                hh(10).FaceColor = [145/255 0 211/255;]; % T3.
                hh(11).FaceColor = [1 0 0]; % T4.
                hh(12).FaceColor = [0 0 0]; % T5.
                P = plot(time,D/1000,'k'); P.LineWidth = 3;
                if UseBattery == 1
                    P = plot(time,(D-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0))/1000,'c');
                    P.LineWidth = 3;
                end
                P = plot(time,(D+Export-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0))/1000,'m');
                P.LineWidth = 1.5;
                if UseBattery == 1
                    if Format_Long == 0
                        legend('Biomass Power','Solar Power','Wind Power','Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)',...
                            'Motores Batlle (Fossil)','PTA (Fossil)','PTB (Fossil)','CTR (Fossil)','Failure (Fossil)','Battery','Demand','Demand + Battery','Demand + Battery + Exp','location','eastoutside');
                    else
                        legend('Biomass Power','Solar Power','Wind Power','Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)',...
                            'Motores Batlle (Fossil)','PTA (Fossil)','PTB (Fossil)','CTR (Fossil)','Failure (Fossil)','Battery','Demand','Demand + Battery','Demand + Battery + Exp','location','eastoutside');
                    end
                else
                    legend('Biomass Power','Solar Power','Wind Power','Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)',...
                        'Motores Batlle (Fossil)','PTA (Fossil)','PTB (Fossil)','CTR (Fossil)','Failure (Fossil)','Demand','Demand + Exp','location','eastoutside');
                end
                xlabel('Time');
                title('Simulated Power Balance')
                ylabel('Power (MW)');
                xlim([0 max(time)]);
                set(findall(gcf,'-property','FontSize'),'FontSize',22);
                set(gca,'Box','on');

                if WhatToDo == 4
                    set(0,'CurrentFigure',101); clf(101); set(gcf,'Position',[10 10 900 420],'name','Power Balance');
                    grid on; grid minor; hold on;
                    hh = area(time,[BioMasPower;SolarPower;WindPower;Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                        Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP]'/1000);
                    hh(1).FaceColor = [1 .5 0]; % Biomasa.
                    hh(2).FaceColor = [1 1 .5]; % Solar.
                    hh(3).FaceColor = [0 1 0]; % Wind.
                    hh(4).FaceColor = [.5 .5 1]; % Bonete.
                    hh(5).FaceColor = [.5 .8 1]; % Baygorria.
                    hh(6).FaceColor = [.2 .2 1]; % Palmar.
                    hh(7).FaceColor = [.8 .8 1]; % SG.
                    hh(8).FaceColor = [1 .6 .6]; % T1.
                    hh(9).FaceColor = [1 .3 .5]; % T2.
                    hh(10).FaceColor = [145/255 0 211/255;]; % T3.
                    hh(11).FaceColor = [1 0 0]; % T4.
                    hh(12).FaceColor = [0 0 0]; % T5.
                    P = plot(time,D/1000,'k'); P.LineWidth = 1;
                    P = plot(time,(D+Export)/1000,'m');
                    P.LineWidth = 1;
                    legend('Biomass Power','Solar Power','Wind Power',Dams{4},Dams{1:3},Termicas{1:5},'Demand','Demand + Exp','location','eastoutside');
                    xlabel('Time');
                    title('Power Balance')
                    ylabel('Power (MW)');
                    xlim([0 max(time)]);
                    set(gca,'Box','on');
                end

                set(0,'CurrentFigure',102); clf(102); set(gcf,'name','Costs Comparison');
                hold on;
                plot(time,Optimal_Cost/NormCost);
                if WhatToDo == 7 || WhatToDo == 4 || WhatToDo == 3 || WhatToDo == 2
                    plot(time,ones(1,length(time))*Final_Cost);
                else
                	plot(time,ones(1,length(time))*Final_Cost/NormCost);
                end
                grid on; grid minor;
                legend('Optimal Path Cost','HJB Final Cost','location','southeast');
                xlabel('Time'); ylabel('Hundred Thousand Dollars');
                title('Costs Comparison');
                xlim([0 max(time)]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',103); clf(103); set(gcf,'name','Controls over the Optimal Path');
                hold on; grid on; grid minor;
                title('Controls over the Optimal Path');
                xlabel('Time');
                ylabel('Control');
                for q = 1:6
                    P = plot(time,Controls_OP(:,q));
                    P.LineWidth = 1;
                end
                legend('Turbine Bonete','Spillage Bonete','Turbine Baygorria','Spillage Baygorria',...
                    'Turbine Palmar','Turbine Salto Grande','location','southeast');
                xlim([0 max(time)]); ylim([0 1]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',104); clf(104); set(gcf,'name','Controls over the Optimal Path');
                hold on; grid on; grid minor;
                title('Controls over the Optimal Path');
                xlabel('Time');
                ylabel('Control');
                for q = 8:12
                    P = plot(time,Controls_OP(:,q));
                    P.LineWidth = 1;
                end
                legend(Termicas{1:4});
                xlim([0 max(time)]); ylim([0 1]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',105); clf(105); set(gcf,'name','Controls over the Optimal Path');
                P = plot(time,Controls_OP(:,7)');
                P.LineWidth = 1; hold on;
                title('Controls over the Optimal Path');
                xlabel('Time'); ylabel('Control');
                plot(time,zeros(1,length(time)),'r');
                legend('Battery');grid on; grid minor;
                xlim([0 max(time)]); ylim([-0.35 1]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',106); clf(106); set(gcf,'name','Volume Bonete');
                plot(time,Vol_OP_1); hold on;
                plot(time,Vol_OP_Real_1);
                ylim([Vmin1 Vmax1]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Bonete');
                legend('Approximated Volume Optimal Path','Real Volume Optimal Path');
                xlim([0 max(time)]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',107); clf(107); set(gcf,'name','Volume Palmar');
                plot(time,Vol_OP_3); hold on;
                plot(time,Vol_OP_Real_3);
                ylim([Vmin3 Vmax3]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Palmar');
                legend('Approximated Volume Optimal Path','Real Volume Optimal Path');
                xlim([0 max(time)]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',108); clf(108); set(gcf,'name','Volume Salto Grande');
                plot(time,Vol_OP_4); hold on;
                plot(time,Vol_OP_Real_4);
                ylim([Vmin4 Vmax4]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Salto Grande');
                legend('Approximated Volume Optimal Path','Real Volume Optimal Path');
                xlim([0 max(time)]);
                set(gca,'Box','on');

                set(0,'CurrentFigure',109); clf(109); set(gcf,'name','Battery Capacity');
                plot(time,Battery_OP); hold on;
                plot(time,Battery_Real);
                ylim([NAMin NAMax]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Capacity');
                title('Battery Capacity');
                legend('Approximated Battery Capacity Path','Real Battery Capacity Path');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',127); clf(127); set(gcf,'name','Battery Capacity');
                P = plot(time,Battery_Real);
                P.LineWidth = 1;
                ylim([NAMin NAMax]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Capacity');
                title('Battery Capacity');
                legend('Battery Capacity over Optimal Path');
                xlim([0 max(time)]);
                set(findall(gcf,'-property','FontSize'),'FontSize',14);
                saveas(gcf,[pwd '/',['Simulation_',name],'/Extra_',num2str(i)],'epsc');
                set(gca,'Box','on');

                set(0,'CurrentFigure',110); clf(110); set(gcf,'name','Costs in the Optimal Path');
                hold on;
                P = plot(time,CostControls(1,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(3,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(5,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(6,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(7,:)*MWhTokW); P.LineWidth = 1;
                grid on; grid minor;
                title('Costs in the Optimal Path');
                legend('Turbine flow Bonete','Turbine flow Baygorria','Turbine flow Palmar','Turbine flow Salto Grande','Battery','location','southeast');
                xlabel('Time'); ylabel('USD/MWh');
                ylim([0.9*min([min(CostControls(1,:)),min(CostControls(3,:)),min(CostControls(5,:)),min(CostControls(6,:))])...
                    1.1*max([max(CostControls(1,:)),max(CostControls(3,:)),max(CostControls(5,:)),max(CostControls(6,:))])]*MWhTokW);
                set(gca,'Box','on');

                set(0,'CurrentFigure',111); clf(111); set(gcf,'name','Costs in the Optimal Path');
                hold on;
                P = plot(time,CostControls(7,:)*MWhTokW);
                P.LineWidth = 1;
                grid on; grid minor;
                title('Costs in the Optimal Path');
                legend('Battery'); xlabel('Time'); ylabel('USD/MWh');
                set(gca,'Box','on');

                set(0,'CurrentFigure',112); clf(112); set(gcf,'name','Costs in the Optimal Path');
                hold on;
                for cont = 8:11
                    P = plot(time,CostControls(cont,:)*MWhTokW);
                    P.LineWidth = 1;
                end
                P = plot(time,CostControls(7,:)*MWhTokW);
                P.LineWidth = 1;
                grid on; grid minor;
                title('Costs in the Optimal Path');
                legend(Termicas{1:4},'Battery'); xlabel('Time'); ylabel('USD/MWh');
                ylim([0.9*min([min(CostControls(8,:)),min(CostControls(9,:)),min(CostControls(10,:)),min(CostControls(11,:))])...
                    1.1*max([max(CostControls(8,:)),max(CostControls(9,:)),max(CostControls(10,:)),max(CostControls(11,:))])]*MWhTokW);
                set(gca,'Box','on');

                set(0,'CurrentFigure',113); clf(113); set(gcf,'name','Costs in the Optimal Path');
                hold on;
                P = plot(time,CostControls(2,:)); P.LineWidth = 1;
                P = plot(time,CostControls(4,:)); P.LineWidth = 1;
                grid on; grid minor;
                title('Costs in the Optimal Path');
                legend('Bonete Spillage','Baygorria Spillage');
                xlabel('Time'); ylabel('USD/hm^3');
                set(gca,'Box','on');

                set(0,'CurrentFigure',114); clf(114); set(gcf,'name','Instant Derivatives Values');
                hold on;
                plot(time,uV1plot/NormCost); plot(time,uV3plot/NormCost);
                plot(time,uV4plot/NormCost); plot(time,uAplot/NormCost);
                grid on; grid minor;
                Leg2 = legend('$u_{V1}$','$u_{V3}$','$u_{V4}$','$u_{A}$','location','southeast');
                set(Leg2,'interpreter', 'latex');
                title('Instant Derivatives Values');
                xlabel('Time');
                ylabel('Derivatives');
                set(gca,'Box','on');

                set(0,'CurrentFigure',115); clf(115); set(gcf,'name','dt/CFL Condition');
                hold on;
                plot(time,CFL_Plot(:,1));
                grid on; grid minor;
                title('CFL Condition');
                ylabel('$\Delta t$/Condition','Interpreter','latex')
                xlabel('Time');
                set(gca,'Box','on');

                set(0,'CurrentFigure',116); clf(116); set(gcf,'name','Energy Distribution');
                X = [sum(BioMasPower),sum(SolarPower),sum(WindPower),sum(Power_H1_OP),sum(Power_H2_OP),...
                    sum(Power_H3_OP),sum(Power_H4_OP),sum(Power_F1_OP),sum(Power_F2_OP),sum(Power_F3_OP),...
                    sum(Power_F4_OP),sum(Power_F5_OP),sum(Power_Battery_OP.*(Power_Battery_OP > 0)),];
                txt = {'Biomass: ','Solar Power: ','Wind Power: ',[Dams{1},': '],[Dams{2},': '],...
                    [Dams{3},': '],[Dams{4},': '],[Termicas{1},': '],[Termicas{2},': '],...
                    [Termicas{3},': '],[Termicas{4},': '],[Termicas{5},': '],'Battery: '};
                colors = [1 .5 0; % Biomasa.
                1 1 .5; % Solar.
                0 1 0; % Wind.
                .5 .5 1;
                .5 .8 1;
                .2 .2 1;
                .8 .8 1;
                1 .6 .6;
                1 .3 .5;
                145/255 0 211/255;
                1 0 0;
                0 0 0;
                1 1 0];
                for i = length(X):-1:1
                    if X(i) <= 0
                        X(i) = [];
                        txt(i) = [];
                        colors(i,:) = [];
                    end
                end
                p1 = pie(X);
                title('Energy Distribution');
                pText = findobj(p1,'Type','text');
                percentValues = get(pText,'String');
                ax = gca;
                delete(ax.Children([1:2:end]));
                combinedtxt = strcat(txt,percentValues');
                legend(combinedtxt,'Location','eastoutside');
                colormap(colors);
                set(gca,'Box','on');

                set(0,'CurrentFigure',117); clf(117); set(gcf,'name','Costs Distribution');
                X = [sum(Cost(1,:)),sum(Cost(2,:)),sum(Cost(3,:)),sum(Cost(4,:)),...
                    sum(Cost(5,:)),sum(Cost(6,:)),sum(Cost(7,:)),sum(Cost(8,:)),sum(Cost(10,:))];
                txt = {[Dams{1},': '],[Dams{2},': '],[Dams{3},': '],[Dams{4},': '],...
                    [Termicas{1},': '],[Termicas{2},': '],[Termicas{3},': '],[Termicas{4},': '],[Termicas{5},': ']};
                colors = [.5 .5 1;
                .5 .8 1;
                .2 .2 1;
                .8 .8 1;
                1 .6 .6;
                1 .3 .5;
                145/255 0 211/255;
                1 0 0
                0 0 0];
                for i = length(X):-1:1
                    if X(i) <= 0
                        X(i) = [];
                        txt(i) = [];
                        colors(i,:) = [];
                    end
                end
                p2 = pie(X);
                title('Controllable Costs Distribution');
                pText = findobj(p2,'Type','text');
                percentValues = get(pText,'String');
                ax = gca;
                delete(ax.Children([1:2:end]));
                combinedtxt = strcat(txt,percentValues');
                legend(combinedtxt,'Location','eastoutside');
                colormap(colors);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',123); clf(123); set(gcf,'name','Costs Histogram');
                X = [sum(Cost(1,:)),sum(Cost(2,:)),sum(Cost(3,:)),sum(Cost(4,:)),...
                    sum(Cost(5,:)),sum(Cost(6,:)),sum(Cost(7,:)),sum(Cost(8,:)),sum(Cost(10,:))]/NormCost;
                txt = {[Dams{1},': ',num2str(round(X(1),3))],[Dams{2},': ',num2str(round(X(2),3))],[Dams{3},': ',num2str(round(X(3),3))],...
                    [Dams{4},': ',num2str(round(X(4),3))],[Termicas{1},': ',num2str(round(X(5),3))],[Termicas{2},': ',num2str(round(X(6),3))],...
                    [Termicas{3},': ',num2str(round(X(7),3))],[Termicas{4},': ',num2str(round(X(8),3))],[Termicas{5},': ',num2str(round(X(9),3))]};
                colors = [.5 .5 1;
                .5 .8 1;
                .2 .2 1;
                .8 .8 1;
                1 .6 .6;
                1 .3 .5;
                145/255 0 211/255;
                1 0 0
                0 0 0];
                p2 = bar(diag(X),'stacked');
                for i = 1:length(X)
                    p2(i).FaceColor = colors(i,:);
                end
                grid minor;
                title('Controllable Costs Distribution');
                legend(txt,'Location','eastoutside');
                ylabel('Hundred Thousand Dollars');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',126); clf(126); set(gcf,'name','Costs Histogram');
                X = [sum(CostR(1,:)),sum(CostR(2,:)),sum(CostR(3,:)),sum(CostR(4,:)),...
                    sum(Cost(5,:)),sum(Cost(6,:)),sum(Cost(7,:)),sum(Cost(8,:)),sum(Cost(10,:))]/NormCost;
                txt = {[Dams{1},': ',num2str(round(X(1),3))],[Dams{2},': ',num2str(round(X(2),3))],[Dams{3},': ',num2str(round(X(3),3))],...
                    [Dams{4},': ',num2str(round(X(4),3))],[Termicas{1},': ',num2str(round(X(5),3))],[Termicas{2},': ',num2str(round(X(6),3))],...
                    [Termicas{3},': ',num2str(round(X(7),3))],[Termicas{4},': ',num2str(round(X(8),3))],[Termicas{5},': ',num2str(round(X(9),3))]};
                colors = [.5 .5 1;
                .5 .8 1;
                .2 .2 1;
                .8 .8 1;
                1 .6 .6;
                1 .3 .5;
                145/255 0 211/255;
                1 0 0
                0 0 0];
                p2 = bar(diag(X),'stacked');
                for i = 1:length(X)
                    p2(i).FaceColor = colors(i,:);
                end
                grid minor;
                title('Costs Distribution');
                legend(txt,'Location','eastoutside');
                ylabel('Hundred Thousand Dollars');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',124); clf(124); set(gcf,'Position',[10 10 900 420],'name','Cost Optimal Path');
                grid minor; hold on;
                hh = area(time,[CostAcum1;CostAcum2;CostAcum3;CostAcum4;CostAcum5;...
                    CostAcum6;CostAcum7;CostAcum8;CostAcum9]'/NormCost);
                hh(1).FaceColor = [.5 .5 1]; % Bonete.
                hh(2).FaceColor = [.5 .8 1]; % Baygorria.
                hh(3).FaceColor = [.2 .2 1]; % Palmar.
                hh(4).FaceColor = [.8 .8 1]; % SG.
                hh(5).FaceColor = [1 .6 .6]; % T1.
                hh(6).FaceColor = [1 .3 .5]; % T2.
                hh(7).FaceColor = [145/255 0 211/255;]; % T3.
                hh(8).FaceColor = [1 0 0]; % T4.
                hh(9).FaceColor = [0 0 0]; % T5.
                txt = [Dams(:)',Termicas(:)'];
                legend(txt,'location','eastoutside');
                xlabel('Time');
                title('Controllable Accumulated Cost')
                ylabel('Hundred Thousand Dollars');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',125); clf(125); set(gcf,'Position',[10 10 900 420],'name','Cost Optimal Path');
                grid minor; hold on;
                hh = area(time,[CostAcum1R;CostAcum2R;CostAcum3R;CostAcum4R;CostAcum5R;...
                    CostAcum6R;CostAcum7R;CostAcum8R;CostAcum9R]'/NormCost);
                hh(1).FaceColor = [.5 .5 1]; % Bonete.
                hh(2).FaceColor = [.5 .8 1]; % Baygorria.
                hh(3).FaceColor = [.2 .2 1]; % Palmar.
                hh(4).FaceColor = [.8 .8 1]; % SG.
                hh(5).FaceColor = [1 .6 .6]; % T1.
                hh(6).FaceColor = [1 .3 .5]; % T2.
                hh(7).FaceColor = [145/255 0 211/255;]; % T3.
                hh(8).FaceColor = [1 0 0]; % T4.
                hh(9).FaceColor = [0 0 0]; % T5.
                P = plot(time,Optimal_Cost/NormCost); P.LineWidth = 3;
                if WhatToDo == 7 || WhatToDo == 4 || WhatToDo == 3
                    P = plot(time,ones(1,length(time))*Final_Cost);
                else
                	P = plot(time,ones(1,length(time))*Final_Cost/NormCost);
                end
                P.LineWidth = 3;
                legend(txt{:},'Optimal Path Cost','HJB Final Cost','location','eastoutside');
                xlabel('Time');
                title('Accumulated Cost')
                ylabel('Hundred Thousand Dollars');
                xlim([0 max(time)]);
                set(gca,'Box','on');

                if WhatToDo ~= 4 && WhatToDo ~= 6 && WhatToDo ~= 7
                    set(0,'CurrentFigure',118); clf(118); set(gcf,'name','Spot Price');
                    plot(time,Spot);grid minor;
                    legend('Spot Price');
                    title('Spot Price');
                    xlabel('Time'); ylabel('USD/MWh');
                    set(gca,'Box','on');
                end

                set(0,'CurrentFigure',119); clf(119); set(gcf,'name','Costs over Optimal Path');
                hold on;
                P = plot(time,CostControlsInstant(1,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControlsInstant(3,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControlsInstant(5,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControlsInstant(6,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControlsInstant(7,:)*MWhTokW); P.LineWidth = 1;
                if WhatToDo ~= 4 && WhatToDo ~= 6 && WhatToDo ~= 7
                P = plot(time,Spot); P.LineWidth = 1;
                end
                grid on; grid minor;
                title('Costs over Optimal Path');
                legend('Turbine flow Bonete','Turbine flow Baygorria','Turbine flow Palmar','Turbine flow Salto Grande','Battery','Spot Price','location','southeast');
                xlabel('Time'); ylabel('USD/MWh');
                ylim([0.9*min([min(CostControls(1,:)),min(CostControls(3,:)),min(CostControls(5,:)),min(CostControls(6,:))])...
                     1.1*max([max(CostControls(1,:)),max(CostControls(3,:)),max(CostControls(5,:)),max(CostControls(6,:))])]*MWhTokW);
                 set(gca,'Box','on');
                
                if Format_Long == 0
                    set(0,'CurrentFigure',120); clf(120); set(gcf,'Position',[10 10 900 420],'name','Power Balance Controllable');
                else
                	set(0,'CurrentFigure',120); clf(120);% set(gcf,'Position',[10 10 750 750],'name','Power Balance Controllable');
                end
                grid on; grid minor; hold on;
                hh = area(time,[Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                    Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP;Power_Battery_OP.*(Power_Battery_OP > 0)]'/1000);
                hh(1).FaceColor = [.5 .5 1]; % Bonete.
                hh(2).FaceColor = [.5 .8 1]; % Baygorria.
                hh(3).FaceColor = [.2 .2 1]; % Palmar.
                hh(4).FaceColor = [.8 .8 1]; % SG.
                hh(5).FaceColor = [1 .6 .6]; % T1.
                hh(6).FaceColor = [1 .3 .5]; % T2.
                hh(7).FaceColor = [145/255 0 211/255;]; % T3.
                hh(8).FaceColor = [1 0 0]; % T4.
                hh(9).FaceColor = [0 0 0]; % T5.
                hh(10).FaceColor = [1 1 0]; % Battery.
                if Format_Long == 0
                    legend(Dams{4},Dams{1:3},Termicas{1:5},'Battery','location','eastoutside');
                else
                    legend(Dams{4},Dams{1:3},Termicas{1:5},'Battery','location','SouthEast');
                end
                xlabel('Time');
                title('Power Balance Controllable')
                ylabel('Power (MW)');
                xlim([0 max(time)]);
                set(findall(gcf,'-property','FontSize'),'FontSize',14);
                set(gca,'Box','on');

                set(0,'CurrentFigure',121); clf(121); set(gcf,'name','Costs Distribution');
                X = [sum(Cost(1,:))+sum(Cost(2,:))+sum(Cost(3,:))+sum(Cost(4,:)),...
                    sum(Cost(5,:))+sum(Cost(6,:))+sum(Cost(7,:))+sum(Cost(8,:))+sum(Cost(10,:))];
                txt = {'Hydropower: ','Fossil fuel power: '};
                colors = [30/255 144/255 1;
                1 69/255 0];
                for i = length(X):-1:1
                    if X(i) <= 0
                        X(i) = [];
                        txt(i) = [];
                        colors(i,:) = [];
                    end
                end
                p2 = pie(X);
                title('Costs Distribution');
                pText = findobj(p2,'Type','text');
                percentValues = get(pText,'String');
                ax = gca;
                delete(ax.Children([1:2:end]));
                combinedtxt = strcat(txt,percentValues');
                legend(combinedtxt,'Location','southoutside');
                colormap(colors);
                set(gca,'Box','on');

                set(0,'CurrentFigure',122); clf(122); set(gcf,'name','Virtual controls over the Optimal Path');
                hold on; grid on; grid minor;
                title('Controls over the Optimal Path');
                xlabel('Time');
                ylabel('Control');
                plot(time,Controls_OP(:,13));
                plot(time,Controls_OP(:,14));
                legend('Virtual Control 1','Virtual Control 2');
                xlim([0 max(time)]); ylim([0 1]);
                set(gca,'Box','on');

                if SavePlots == 1 && WhatToDo ~= 4
                    for i = Plot_12:Plot_22
                        set(0,'CurrentFigure',i);
                        saveas(gcf,[pwd '/',['Simulation_',name],'/',num2str(i)],'epsc');
                    end
                end
                if WhatToDo == 4
                    sets = [102,103,104,110,112,116,117];
                    for m = 1:length(sets)
                        i = sets(m);
                        set(0,'CurrentFigure',i);
                        saveas(gcf,[pwd '/',['Simulation_',name],'/R_',num2str(i)],'epsc');
                    end
                end

            end

            Final_Cost = Optimal_Cost(end)/NormCost;
            if GradNum == 1000 
                save([pwd '/','Simulation_Controls_',name,'_Day_',num2str(Day),'.mat'],'Controls_OP');
            end
            toc

            if WhatToDo == 6

                Delta = (length(time)-Z21-1)/(3*2^GradNum);

                for j = 0:3*2^GradNum-1
                    Final_Cost(j+1) = sum(WaterVirtualBay(Z21+1+Delta*j:Z21+Delta*(1+j))-WaterBonete(1+Delta*j:Delta*(1+j)));
                    DiffWater(Z21+1+Delta*j:Z21+Delta*(1+j)) = Final_Cost(j+1);
                end
                DiffWater(end+1) = DiffWater(end);
                Final_Cost = Final_Cost/1e6;

                figure(1000);
                set(0,'CurrentFigure',1000); clf(1000); set(gcf,'name','Gradient');
                plot(time,DiffWater/1e6);
                grid minor;
                title('Difference in Water');
                xlabel('Time');
                ylabel('hm^3');
                xlim([0 max(time)]);
                set(gca,'Box','on');

                figure(1001); set(gcf,'name','All Gradient');
                plot(time,DiffWater/1e6);
                hold on; 
                title('All Differences in Water');
                xlabel('Time');
                ylabel('hm^3');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                grid minor;

                figure(1002);
                set(0,'CurrentFigure',1002); clf(1002); set(gcf,'name','Lagrangian Multipliers');
                lamTime = linspace(0,1+t32,length([Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]));
                plot(lamTime,[Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]*1e6);
                grid minor;
                title('Approximated Lagrangian Multipliers');
                ylabel('USD/hm$^3$','Interpreter','latex');
                xlabel('Time');
                Leg = legend('$\hat{\lambda}_{2,1}(t)$');
                set(Leg, 'interpreter', 'latex');
                xlim([0 max(lamTime)]);

            end
        
        end
        
    end
    
    %% End.

end