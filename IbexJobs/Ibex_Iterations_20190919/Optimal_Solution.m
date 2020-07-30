function [Final_Cost] = Optimal_Solution(name,WhatToDo,WarmStart,UseBattery,...
    ComputePlots,NewPlots,SavePlots,Expan,Expan_Time,UseLamb21,UseLamb32,...
    Lambda21,Lambda32,Day,GradNum)

% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email: renzo.caballerorosas@kaust.edu.sa caballerorenzo@hotmail.com
% Website: None.
% June 2019; Last revision: 10/10/2019.

% name is a string, if name = 'example' then the simulations will save in ./Simulations/Simulation_example.
% Day is a string, if Day = '20190101' then we load Day_20190101.

% WhatToDo can have the next values:
% 1 - We solve the HJB equation.
% 2 - We compute the optimal controls.
% 3 - We do 1 and after 2.
% 4 - 
% 5 - 
% 6 - We compute the optimal controls but letting the virtual controls to
% have any value, and also removing the H1 penalization.
% 7 - We load and use the historical controls.
%
% When GradNum = 1000, we save the controls of the optimal patha after we
% compute them.
% The controls are [turBon, spillBon, turBay, spillBay, turPal,
% turSA, Battery, Mot, PTA, PTB, CTR, Failure, V2, V3].

%     Lambda21 = Lambda21/1e6; % We pass from USD/hm^3 to USD/m^3.
%     Lambda32 = Lambda32/1e6; % We pass from USD/hm^3 to USD/m^3.
    % GradNum is in {0,1,...,n-2} where n is Expan_Time.
    Derivatives = [1,1,1,1,1]; % For debugging.
    artificialWaterValue = 0; % For debugging.
    artificialFFS = 0; % For debugging.
    artificialVirtualControls = 0; % For debugging.
    postWithZeroPartials = 0; % For debugging.
    
    % Remmeber to vectorize everything, i.e., V(v1,v2,...,vn) = V(vec) with
    % vec = [v1,v2,...,vn].
    
    if 0 == exist([pwd '/Simulations/Simulation_',name],'dir')
        mkdir([pwd '/Simulations/Simulation_',name]);
    end

    load('MaxFlowCoeff.mat');
    load([pwd '/Historical/dailyData/Day_',Day,'.mat']);
    
    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',1e4,...
        'FunctionTolerance',1e-10,'MaxIterations',1e3,'MaxFunEvals',1e7,'Algorithm','sqp',...
        'ConstraintTolerance',1e-5,'StepTolerance',1e-15);

    Plot_11 = 1;
    Plot_21 = 15;
    Plot_12 = 101;
    Plot_22 = 133;
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
    Termicas = {'Motores Batlle','PTA','PTB','CTR','Failure'};
    Dams = {'Bonete','Baygorria','Palmar','Salto Grande'};
    MaxTermicas = [Matrix{2}(4),Matrix{2}(6),Matrix{2}(8),Matrix{2}(2),1000]*1000; % kW.
    InfTermicas = [0,0,0,0,0];
    
    CostTermicas = [Matrix{2}(3),Matrix{2}(5),Matrix{2}(7),Matrix{2}(1),600]; % USD/MWh.
    if artificialFFS
%         CostTermicas = CostTermicas/1000;
        CostTermicas(1) = 10;
        CostTermicas(2) = 10;
        CostTermicas(3) = 10;
        CostTermicas(4) = 10;
        CostTermicas(5) = 10;
    end
    CostTermicasNorm = CostTermicas/MWhTokW; % USD/kW.

    NT = 2^Expan_Time; % Discretizations in time.
    NV1 = 2^Expan; % Discretizations of Bonete.
    NV2 = 2^Expan; % Discretizations of Baygorria.
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
    T = length(time);
    
    % ==================== Parameters (Demand, Water and Fuel) ====================>>
       
    H1 = @(V1) (-3.77)*(V1).^2 + (15.7)*(V1) + 69.6; % Height of the water as a function of the normal volume (Bonete).
    H2 = @(V2) (-1.40)*(V2).^2 + (8.89)*(V2) + 47.5;
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
    
    dH1 = @(V1,V2) H1(V1) - H2(V2) - dBB;
    dH2 = @(V2,V3) H2(V2) - H3(V3) - dBP;
    dH3 = @(V3,t) H3(V3) - h03(t);
    dH4 = @(V4,t) H4(V4) - h04(t);
    
    Co1 = MaxFlowCoeff{1};
    Co2 = MaxFlowCoeff{2};
    Co3 = MaxFlowCoeff{3};
    
    % THIS IS A MODIFICATION: To make fixed the dependence between.
    % Otherwise, the simulation is complicated.
    VolMax2 = 0.678e9; % In m^3.
    Vini2 = Matrix{5}(1)/VolMax2*1000000;
    dH1 = @(V1,V2) H1(V1) - H2(Vini2) - dBB;
    VolMax3 = 3.53e9; % In m^3.
    Vini3 = Matrix{6}(1)/VolMax3*1000000;
    dH2 = @(V2,V3) H2(V2) - H3(Vini3) - dBP;

    FlowMax1 = @(V1,V2) Co1(1)*dH1(V1,V2).^2 + Co1(2)*dH1(V1,V2) + Co1(3); % In m^3/s.
    FlowMax2 = @(V2,V3) (dH2(V2,V3) <= Co2{3}).*(Co2{1}(1)*dH2(V2,V3) + Co2{1}(2)) + ...
        (dH2(V2,V3) > Co2{3}).*(Co2{2}(1)*dH2(V2,V3) + Co2{2}(2)); % In m^3/s.
    FlowMax3 = @(V3,t) (dH3(V3,t) <= Co3{3}).*(Co3{1}(1)*dH3(V3,t) + Co3{1}(2)) + ...
        (dH3(V3,t) > Co3{3}).*(Co3{2}(1).*dH3(V3,t).^2 + Co3{2}(2).*dH3(V3,t) + Co3{2}(3)); % In m^3/s.
    FlowMax4 = 4410; % In m^3/s.
    
    Vols = 0:0.01:1;
    RealMaxFlow1 = max(FlowMax1(Vols,Vols)); % In m^3/s.
    RealMaxFlow2 = max(FlowMax2(Vols,Vols)); % In m^3/s.
    RealMaxFlow3 = max(FlowMax3(Vols,find(h03_Vec==max(h03_Vec), 1 ))); % In m^3/s.
    % Just above, we added a min(...) because sometimes
    % h03_Vec == max(h03_Vec) has multiple solutions and the function only
    % admits an integer.
    RealMaxFlow4 = FlowMax4; % In m^3/s.
    
    MaxTot1 = RealMaxFlow1*2; % In m^3/s.
    MaxTot2 = RealMaxFlow2*2; % In m^3/s.
    MaxTot3 = RealMaxFlow3*2; % In m^3/s.
    MaxTot4 = RealMaxFlow4*2; % In m^3/s.
    
    MaxSpill1 = @(V1,V2) MaxTot1 - FlowMax1(V1,V2); % In m^3/s.
    MaxSpill2 = @(V2,V3) MaxTot2 - FlowMax2(V2,V3); % In m^3/s.
    MaxSpill3 = @(V3,t) MaxTot3 - FlowMax3(V3,t); % In m^3/s.
    MaxSpill4 = MaxTot4 - RealMaxFlow4; % In m^3/s.

    MaxV2 = MaxTot1; % In m^3/s.
    MaxV3 = MaxTot2; % In m^3/s.

    VolMax1 = 10.7e9; % In m^3.
    VolMax2 = 0.678e9; % In m^3.
    VolMax3 = 3.53e9; % In m^3.
    VolMax4 = 5.18e9; % In m^3.

    VolMin1 = 1.85e9; % In m^3. From tables at 70 m.
    VolMin2 = 0.47e9;
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
    IT4 = @(t) Inflow4/2; % In m^3/s (we divide because Uruguay and Argentina share the water).

    f1 = @(t,V1,V2,T,S) (TMax/VolMax1)*(IT1(t)-FlowMax1(V1,V2)*T-MaxSpill1(V1,V2)*S);
    f2 = @(t,V2,V3,T,S,V) (TMax/VolMax2)*(IT2(t)-FlowMax2(V2,V3)*T-MaxSpill2(V2,V3)*S+MaxV2*V);
    f3 = @(t,V3,T,S,V) (TMax/VolMax3)*(IT3(t)-FlowMax3(V3,t)*T-MaxSpill3(V3,t)*S+MaxV3*V);
    f4 = @(t,V4,T,S) (TMax/VolMax4)*(IT4(t)-FlowMax4*T-MaxSpill4*S);

    DMax = 2.5e6; % In kW. The maximum possible value of the demand.

    Vini1 = Matrix{4}(1)/VolMax1*1000000;
    Vini2 = Matrix{5}(1)/VolMax2*1000000;
    Vini3 = Matrix{6}(1)/VolMax3*1000000;
    Vini4 = VolSG(Matrix{7}(2));

    t21 = 1/4; % Delay between Bonete and Baygorria.
    Z21 = floor(length(time)*t21-1)+1; % Discrete delay between Bonete and Baygorria.
    t32 = 10/24; % Delay between Baygorria and Palmar.
    Z32 = floor(length(time)*t32-1)+1; % Discrete delay between Baygorria and Palmar.
    
    Phi_V2 = Matrix{9}(:,1)/MaxV2; % We divide to normalize it (it becomes between 0 and 1).
    Phi_V3 = Matrix{9}(:,2)/MaxV3; % We divide to normalize it.
    Phi_V2 = interp1(linspace(0,1,length(Phi_V2)),Phi_V2,linspace(0,1,Z21-1));
    Phi_V3 = interp1(linspace(0,1,length(Phi_V3)),Phi_V3,linspace(0,1,Z32-1));

    if UseLamb21 == 0 || (UseLamb21 == 1 && length(Lambda21) == 1)
        Lambda21 = 0; % In USD/m^3.
    end
    if UseLamb32 == 0 || (UseLamb32 == 1 && length(Lambda32) == 1)
        Lambda32 = 0; % In USD/m^3.
    end
    lenLambda21 = length(Lambda21);
    lenLambda32 = length(Lambda32);
    D21         = floor((T-Z21)/lenLambda21);
    D32         = floor((T-Z32)/lenLambda32);
    auxLambda21 = Lambda21*1e-6; % We pass from USD/hm^3 to USD/m^3.
    auxLambda32 = Lambda32*1e-6; % We pass from USD/hm^3 to USD/m^3.
    Lambda21    = zeros(1,length(time)+max(Z21,Z32));
    Lambda32    = zeros(1,length(time)+max(Z21,Z32));
    for i = 0:lenLambda21-1
        if i == lenLambda21-1
            Lambda21(Z21+i*D21:T-1) = auxLambda21(i+1);
        else
            Lambda21(Z21+i*D21:Z21+(i+1)*D21-1) = auxLambda21(i+1);
        end
    end
    for i = 0:lenLambda32-1
        if i == lenLambda32-1
            Lambda32(Z32+i*D32:T-1) = auxLambda32(i+1);
        else
            Lambda32(Z32+i*D32:Z32+(i+1)*D32-1) = auxLambda32(i+1);
        end
    end
    
    if ComputePlots == 1
        set(0,'CurrentFigure',11); clf(11); set(gcf,'name','Lagrangian Multipliers');
        lamTime = linspace(0,1+t32,length(Lambda21));
        plot(lamTime,Lambda21 * 1e6); % We multiply by 1e6 to pass from USD/m^3 to USD/hm^3.
        hold on; grid on;
        plot(lamTime,Lambda32 * 1e6); % We multiply by 1e6 to pass from USD/m^3 to USD/hm^3.
        title('Approximated Lagrangian Multipliers');
        ylabel('USD/hm$^3$','Interpreter','latex');
        xlabel('Time');
        legend('$\hat{\lambda}_{2,1}(t)$','$\hat{\lambda}_{3,2}(t)$','interpreter','latex');
        xlim([0 max(lamTime)]);
    elseif ComputePlots ~= 0
        error('Choose ComputePlots between 0 or 1.');
    end

    NormCost = 1e5; % The normalization is over 100.000 dollars.
    expFactor = 1.5; % This factor is to extend the domain of the dams.
    
    Vmax1 = min((Vini1*VolMax1+(MaxTot1)*TMax*expFactor),VolMax1)/VolMax1;
    Vmin1 = max((Vini1*VolMax1-(MaxTot1)*TMax*expFactor),VolMin1)/VolMax1;
    Vmax2 = min((Vini2*VolMax2+(MaxTot2)*TMax*expFactor),VolMax2)/VolMax2;
    Vmin2 = max((Vini2*VolMax2-(MaxTot2)*TMax*expFactor),VolMin2)/VolMax2;
    Vmax3 = min((Vini3*VolMax3+(RealMaxFlow3)*TMax*expFactor),VolMax3)/VolMax3;
    Vmin3 = max((Vini3*VolMax3-(RealMaxFlow3)*TMax*expFactor),VolMin3)/VolMax3;
    Vmax4 = min((Vini4*VolMax4+(RealMaxFlow4)*TMax*expFactor),VolMax4)/VolMax4;
    Vmin4 = max((Vini4*VolMax4-(RealMaxFlow4)*TMax*expFactor),VolMin4)/VolMax4;
    
    dV1 = (Vmax1-Vmin1)/NV1;
    dV2 = (Vmax2-Vmin2)/NV2;
    dV3 = (Vmax3-Vmin3)/NV3;
    dV4 = (Vmax4-Vmin4)/NV4;

    V1 = Vmin1:dV1:Vmax1; % Discretized volume of Bonete.
    V2 = Vmin2:dV2:Vmax2; % Discretized volume of Baygorria.
    V3 = Vmin3:dV3:Vmax3; % Discretized volume of Palmar.
    V4 = Vmin4:dV4:Vmax4; % Discretized volume of Salto Grande.

%     [vplot1_4,vplot3_4] = meshgrid(V1,V3); % For plotting surfaces.
%     [vplot1_3,vplot4_3] = meshgrid(V1,V4); % For plotting surfaces.
    [vplot3_1,vplot4_1] = meshgrid(V3,V4); % For plotting surfaces.
    
    D = Matrix{1}(:,12)'*1000;%+400e3;
    
    Export = (Matrix{1}(:,84)+Matrix{1}(:,85)+Matrix{1}(:,86))*1000;
    WindPower = Matrix{1}(:,5)'*1000;
    SolarPower = Matrix{1}(:,6)'*1000;
    BioMasPower = Matrix{1}(:,8)'*1000;
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
        legend('Demand','Demand + Exportation');
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
        hold on; grid on; grid minor;
        xlabel('A');
        ylabel('Control');
        title('Battery Control Vs Charge');
        P = plot(A,olPn(A)); P.LineWidth = 1;
        P = plot(A,ulPn(A)); P.LineWidth = 1;
        legend('Maximum Control','Minimum Control');
        pause(0.001);
    end

%     [vplotA_1,vplot1_A] = meshgrid(A,V1); % For plotting surfaces.
    [vplotA_3,vplot3_A] = meshgrid(A,V3); % For plotting surfaces.
    [vplotA_4,vplot4_A] = meshgrid(A,V4); % For plotting surfaces.

    fA = @(C) TMax*K_A*C;
    KA = 0.3/(MWhTokW); % In USD/kW.

    % ==================== Matrices ====================>>

    u = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    CFL = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    CFL_v1 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    CFL_v2 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    CFL_v3 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    CFL_v4 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    CFL_A = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    HAM = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    uV1 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    uV2 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    uV3 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    uV4 = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    uA = zeros(length(V1),length(V2),length(V3),length(V4),length(A));
    ControlsFMC = zeros(length(V1),length(V2),length(V3),length(V4),length(A),14); % T1,S1,T2,S2,T3,T4,A,F1,F2,F3,F4,V2,V3.
    LV1 = length(V1); % Length of V1.
    LV2 = length(V2); % Length of V1.
    LV3 = length(V3); % Length of V3.
    LV4 = length(V4); % Length of V4.
    LA = length(A); % Length of A.
    TL = LV1*LV2*LV3*LV4*LA; % Total number of elements in the matrix.
    Controls_T = ControlsFMC;
    UV = {};
    U = {};

    % ==================== More parameters (Water and Fuel) ====================>>

    c1 = -MaxFlowCoeff{4}{1};
    eta1 = MaxFlowCoeff{4}{2};
    d1 = @(V1,V2,tur,spill) c1*(tur*FlowMax1(V1,V2)+ spill*MaxSpill1(V1,V2));
    % Variation in the high due of the released flow (Bonete).
    d1NonNorm = @(tur,spill) c1*(tur + spill); % In this case, the flows are in m^3/s.

    c2 = -MaxFlowCoeff{4}{3};
    eta2 = MaxFlowCoeff{4}{4};
    d2 = @(V2,V3,tur,spill) c2*(tur*FlowMax2(V2,V3)+spill*MaxSpill2(V2,V3));
    % Variation in the high due of the released flow (Baygorria).
    d2NonNorm = @(tur,spill) c2*(tur + spill); % In this case, the flows are in m^3/s.

    c3 = -MaxFlowCoeff{4}{5};
    eta3 = MaxFlowCoeff{4}{6};
    d3 = @(V3,tur,spill,t) c3*(tur*FlowMax3(V3,t)+spill*MaxSpill3(V3,t));
    % Variation in the high due of the released flow (Palmar).

    c4 = 0.00147;
    eta4 = 8.82*1.16; % The 1.16 is a correction factor, related with the data of real production.
    d4 = @(tur,spill) c4*(tur*FlowMax4+spill*MaxSpill4);
    % Variation in the high due of the released flow (Salto Grande).

    Khe1 = Matrix{3}(3); % In USD/MWh. Cost per energy of Bonete.
    Khe2 = Matrix{3}(4) + 0.1; % In USD/MWh. Cost per energy of Baygorria.
    Khe3 = Matrix{3}(2); % In USD/MWh. Cost per energy of Palmar.
    Khe4 = Matrix{3}(1); % In USD/MWh. Cost per energy of SG.
    
    if artificialWaterValue
%         Khe1 = 0; Khe2 = 0; Khe3 = 0; Khe4 = 0; 
%         Khe1 = 1;
        Khe2 = 0;
        Khe3 = 0;
%         Khe4 = 1; 
    end

    CH1 = Khe1*eta1*(dH1(Vini1,Vini2)-d1(Vini1,Vini2,1,0))/(MWhTokW); % Water's value of Bonete (USD/m^3).
    CH2 = Khe2*eta2*(dH2(Vini2,Vini3)-d2(Vini2,Vini3,1,0))/(MWhTokW); % Water's value of Baygorria (USD/m^3).
    CH3 = Khe3*eta3*(dH3(Vini3,1)-d3(Vini3,1,0,1))/(MWhTokW); % Water's value of Palmar (USD/m^3).
    CH4 = Khe4*eta4*(dH4(Vini4,1)-d4(1,0))/(MWhTokW); % Water's value of Salto Grande (USD/m^3).

    olCT1 = @(uV1,t,V1,V2) FlowMax1(V1,V2)*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCS1 = @(uV1,t,V1,V2) MaxSpill1(V1,V2)*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCT2 = @(uV2,t,V2,V3) FlowMax2(V2,V3)*(CH2 - Lambda32(t+Z32) - uV2/VolMax2);
    olCS2 = @(uV2,t,V2,V3) MaxSpill2(V2,V3)*(CH2 - Lambda32(t+Z32) - uV2/VolMax2);
    olCT3 = @(uV3,t,V3)    FlowMax3(V3,t)*(CH3 - uV3/VolMax3);
    olCT4 = @(uV4)         FlowMax4*(CH4 - uV4/VolMax4);
    olCA  = @(uA)          -K_A*uA;
    olCV2 = @(uV2,t)       MaxV2*(Lambda21(t) + uV2/VolMax2 - CH2);
    olCV3 = @(uV3,t)       MaxV3*(Lambda32(t) + uV3/VolMax3 - CH3);
    H0_1  = @(uV1,t)       IT1(t)*(uV1/VolMax1-CH1);
    H0_2  = @(uV2,t)       IT2(t)*(uV2/VolMax2-CH2);
    H0_3  = @(uV3,t)       IT3(t)*(uV3/VolMax3-CH3);
    H0_4  = @(uV4,t)       IT4(t)*(uV4/VolMax4-CH4);
    
    H0 = @(uV1,uV2,uV3,uV4,t) H0_1(uV1,t)+H0_2(uV2,t)+H0_3(uV3,t)+H0_4(uV4,t);
    
%     k1 = @(t) IT2(t)/MaxV2;
%     k2 = @(V2,V3) FlowMax2(V2,V3)/MaxV2;
%     k3 = @(V2,V3) MaxSpill2(V2,V3)/MaxV2;
%     k = @(t,V2,V3) [k1(t),k2(V2,V3),k3(V2,V3)];

    % Hydraulic power: PH = T^2*K2 + T*K1 + T*S*K3.

    K11 = @(V1,V2) eta1*FlowMax1(V1,V2).*dH1(V1,V2);
    K12 = @(V2,V3) eta2*FlowMax2(V2,V3).*dH2(V2,V3);
    K13 = @(V3,t) eta3*FlowMax3(V3,t).*dH3(V3,t);
    K14 = @(V4,t) eta4*FlowMax4*dH4(V4,t);

    K21 = @(V1,V2) -eta1*FlowMax1(V1,V2).*d1(V1,V2,1,0);
    K22 = @(V2,V3) -eta2*FlowMax2(V2,V3).*d2(V2,V3,1,0);
    K23 = @(V3,t) -eta3*FlowMax3(V3,t).*d3(V3,1,0,t);
    K24 = -eta4*FlowMax4.*d4(1,0);
    
    K31 = @(V1,V2) -eta1*FlowMax1(V1,V2)*d1(V1,V2,0,1);
    K32 = @(V2,V3) -eta2*FlowMax2(V2,V3)*d2(V2,V3,0,1);

    d = @(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t)...
        [olCT1(uV1,t,V1,V2),...
        olCS1(uV1,t,V1,V2),...
        olCT2(uV2,t,V2,V3),...
        olCS2(uV2,t,V2,V3),...
        olCT3(uV3,t,V3),...
        olCT4(uV4),...
        olCA(uA),...
        CostTermicasNorm.*MaxTermicas]';
    dd = diag([0,0,0,0,0,0,KA*1000000*0.000,zeros(1,length(CostTermicasNorm))]);
    
    % We want to minimize min(x^t.d+x^t.dd.x), dd is the quadratic cost. d is
    % the linear cost, and x are the controls.
    
    rho = KA*1000000*0.1*0;
    
    if WhatToDo == 6
        % We remove the penalizations when WhatToDo == 6.
        rho = 0;
    end
    
    % The constrain is x^T.Q.x + x^T.b - c = 0 (or >=).
    b = @(V1,V2,V3,V4,t) [K11(V1,V2),0,K12(V2,V3),0,K13(V3,t),K14(V4,t),PAMax,MaxTermicas]';
    
    Q = @(V1,V2,V3,t) [K21(V1,V2) K31(V1,V2)/2 0 0 0 0 0 0 0 0 0 0;
        K31(V1,V2)/2 0 0 0 0 0 0 0 0 0 0 0;
        0 0 K22(V2,V3) K32(V2,V3)/2 0 0 0 0 0 0 0 0;
        0 0 K32(V2,V3)/2 0 0 0 0 0 0 0 0 0;
        0 0 0 0 K23(V3,t) 0 0 0 0 0 0 0;
        0 0 0 0 0 K24 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0]; % Is fixed.

    c = @(t) D(t)+Export(t)-WindPower(t)-SolarPower(t)-BioMasPower(t);
    
    disp(['KH1 = ',num2str(Khe1),', KH2 = ',num2str(Khe2),', KH3 = ',num2str(Khe3),', KH4 = ',num2str(Khe4),...
        ', KF1 = ',num2str(CostTermicas(1)),', KF2 = ',num2str(CostTermicas(2)),', KF3 = ',...
        num2str(CostTermicas(3)),', KF4 = ',num2str(CostTermicas(4)),', KF5 = ',num2str(CostTermicas(5))]);

    % ==================== Simulation ====================>>
        
    if WhatToDo == 1 || WhatToDo == 3

        tic

        for t = length(time):-1:1

            ControlsFMC_Aux = Controls_T;
            CFL(:) = 0;
            CFL_v1(:) = 0;
            CFL_v2(:) = 0;
            CFL_v3(:) = 0;
            CFL_v4(:) = 0;
            CFL_A(:) = 0;
            
            % ==================== Finite Differences Section ====================>>

            for v1 = 1:length(V1)
                for v2 = 1:length(V3)
                    for v3 = 1:length(V3)
                        for v4 = 1:length(V4)
                            for a = 1:length(A)

                                if (v1~=1) && (v2~=1) && (v3~=1) && (v4~=1) && (a~=1) && (v1~=length(V1)) && (v2~=length(V2)) &&...
                                        (v3~=length(V3)) && (v4~=length(V4)) && (a~=length(A)) && (t<length(time)-1)

                                    if (IT1(t+1)-FlowMax1(V1(v1),V2(v2))*ControlsFMC_Aux(v1,v2,v3,v4,a,1)-MaxSpill1(V1(v1),V2(v2))*ControlsFMC_Aux(v1,v2,v3,v4,a,2) >= 0) && Derivatives(1)
                                        uV1(v1,v2,v3,v4,a) = (u(v1+1,v2,v3,v4,a)-u(v1,v2,v3,v4,a))/(dV1);
                                    elseif Derivatives(1)
                                        uV1(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1-1,v2,v3,v4,a))/(dV1);
                                    end
                                    
                                    if (IT2(t+1)-FlowMax2(V2(v2),V3(v3))*ControlsFMC_Aux(v1,v2,v3,v4,a,3)-MaxSpill2(V2(v2),V3(v3))*ControlsFMC_Aux(v1,v2,v3,v4,a,4)+MaxV2*ControlsFMC_Aux(v1,v2,v3,v4,a,Pc2) >= 0) && Derivatives(2)
                                        uV2(v1,v2,v3,v4,a) = (u(v1,v2+1,v3,v4,a)-u(v1,v2,v3,v4,a))/(dV2);
                                    elseif Derivatives(1)
                                        uV2(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2-1,v3,v4,a))/(dV2);
                                    end

                                    if (IT3(t+1)-FlowMax3(V3(v3),t+1)*ControlsFMC_Aux(v1,v2,v3,v4,a,5)+MaxV3*ControlsFMC_Aux(v1,v2,v3,v4,a,Pc3) >= 0) && Derivatives(3)
                                        uV3(v1,v2,v3,v4,a) = (u(v1,v2,v3+1,v4,a)-u(v1,v2,v3,v4,a))/(dV3);
                                    elseif Derivatives(2)
                                        uV3(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2,v3-1,v4,a))/(dV3);
                                    end

                                    if (IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v2,v3,v4,a,6) >= 0) && Derivatives(4)
                                        uV4(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4+1,a)-u(v1,v2,v3,v4,a))/(dV4);
                                    elseif Derivatives(3)
                                        uV4(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2,v3,v4-1,a))/(dV4);
                                    end

                                    if (-K_A*ControlsFMC_Aux(v1,v2,v3,v4,a,7) >= 0) && Derivatives(5)
                                        uA(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a+1)-u(v1,v2,v3,v4,a))/(dA);
                                    elseif Derivatives(4)
                                        uA(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2,v3,v4,a-1))/(dA);
                                    end

                                elseif t < length(time) - 1

                                    if (v2==1 || v2==length(V2) ||v3==1 || v3==length(V3) || v4==1 || v4==length(V4) || a==1 || a==length(A)) && (v1~=1) && (v1~=length(V1)) && Derivatives(1)
                                        if IT1(t+1)-FlowMax1(V1(v1),V2(v2))*ControlsFMC_Aux(v1,v2,v3,v4,a,1)-MaxSpill1(V1(v1),V2(v2))*ControlsFMC_Aux(v1,v2,v3,v4,a,2) >= 0
                                            uV1(v1,v2,v3,v4,a) = (u(v1+1,v2,v3,v4,a)-u(v1,v2,v3,v4,a))/(dV1);
                                        else
                                            uV1(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1-1,v2,v3,v4,a))/(dV1);
                                        end
                                    end
                                    
                                    if (v1==1 || v1==length(V1) ||v3==1 || v3==length(V3) || v4==1 || v4==length(V4) || a==1 || a==length(A)) && (v2~=1) && (v2~=length(V2)) && Derivatives(2)
                                        if IT2(t+1)-FlowMax2(V2(v2),V3(v3))*ControlsFMC_Aux(v1,v2,v3,v4,a,3)-MaxSpill2(V2(v2),V3(v3))*ControlsFMC_Aux(v1,v2,v3,v4,a,4)+MaxV2*ControlsFMC_Aux(v1,v2,v3,v4,a,Pc2) >= 0
                                            uV2(v1,v2,v3,v4,a) = (u(v1,v2+1,v3,v4,a)-u(v1,v2,v3,v4,a))/(dV2);
                                        elseif Derivatives(1)
                                            uV2(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2-1,v3,v4,a))/(dV2);
                                        end
                                    end

                                    if (v1==1 || v1==length(V1) || v2==1 || v2==length(V2) || v4==1 || v4==length(V4) || a==1 || a==length(A)) && (v3~=1) && (v3~=length(V3)) && Derivatives(3)
                                        if IT3(t+1)-FlowMax3(V3(v3),t+1)*ControlsFMC_Aux(v1,v2,v3,v4,a,5)+MaxV3*ControlsFMC_Aux(v1,v2,v3,v4,a,Pc3) >= 0
                                            uV3(v1,v2,v3,v4,a) = (u(v1,v2,v3+1,v4,a)-u(v1,v2,v3,v4,a))/(dV3);
                                        else
                                            uV3(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2,v3-1,v4,a))/(dV3);
                                        end
                                    end

                                    if (v1==1 || v1==length(V1) || v2==1 || v2==length(V2) || v3==1 || v3==length(V3) || a==1 || a==length(A)) && (v4~=1) && (v4~=length(V4)) && Derivatives(4)
                                        if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v2,v3,v4,a,6) >= 0
                                            uV4(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4+1,a)-u(v1,v2,v3,v4,a))/(dV4);
                                        else
                                            uV4(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2,v3,v4-1,a))/(dV4);
                                        end
                                    end

                                    if (v1==1 || v1==length(V1) || v2==1 || v2==length(V2) || v3==1 || v3==length(V3) || v4==1 || v4==length(V4)) && (a~=1) && (a~=length(A)) && Derivatives(5)
                                        if -K_A*ControlsFMC_Aux(v1,v2,v3,v4,a,7) >= 0
                                            uA(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a+1)-u(v1,v2,v3,v4,a))/(dA);
                                        else
                                            uA(v1,v2,v3,v4,a) = (u(v1,v2,v3,v4,a)-u(v1,v2,v3,v4,a-1))/(dA);
                                        end
                                    end

                                end

                                CFL(v1,v2,v3,v4,a) = 1 / (abs(f1(t,V1(v1),V2(v2),ControlsFMC_Aux(v1,v2,v3,v4,a,1),ControlsFMC_Aux(v1,v2,v3,v4,a,2)))/dV1 + ...
                                    abs(f2(t,V2(v2),V3(v3),ControlsFMC_Aux(v1,v2,v3,v4,a,1),ControlsFMC_Aux(v1,v2,v3,v4,a,2),ControlsFMC_Aux(v1,v2,v3,v4,a,Pc2)))/dV2 + ...
                                    abs(f3(t,V3(v3),ControlsFMC_Aux(v1,v2,v3,v4,a,5),0,ControlsFMC_Aux(v1,v2,v3,v4,a,Pc3)))/dV3 + ...
                                    abs(f4(t,V4(v4),ControlsFMC_Aux(v1,v2,v3,v4,a,6),0))/dV4 + ...
                                    abs(fA(ControlsFMC_Aux(v1,v2,v3,v4,a,7)))/dA);
                                CFL_v1(v1,v2,v3,v4,a) = abs(f1(t,V1(v1),V2(v2),ControlsFMC_Aux(v1,v2,v3,v4,a,1),ControlsFMC_Aux(v1,v2,v3,v4,a,2)))/dV1;
                                CFL_v2(v1,v2,v3,v4,a) = abs(f2(t,V2(v2),V3(v3),ControlsFMC_Aux(v1,v2,v3,v4,a,1),ControlsFMC_Aux(v1,v2,v3,v4,a,2),ControlsFMC_Aux(v1,v2,v3,v4,a,Pc2)))/dV2;
                                CFL_v3(v1,v2,v3,v4,a) = abs(f3(t,V3(v3),ControlsFMC_Aux(v1,v2,v3,v4,a,5),0,ControlsFMC_Aux(v1,v2,v3,v4,a,Pc3)))/dV3;
                                CFL_v4(v1,v2,v3,v4,a) = abs(f4(t,V4(v4),ControlsFMC_Aux(v1,v2,v3,v4,a,6),0))/dV4;
                                CFL_A(v1,v2,v3,v4,a) = abs(fA(ControlsFMC_Aux(v1,v2,v3,v4,a,7)))/dA;

                            end
                        end
                    end
                end
            end

            uV1(1,:,:,:,:) = uV1(2,:,:,:,:);
            uV1(end,:,:,:,:) = uV1(end-1,:,:,:,:);
            uV2(:,1,:,:,:) = uV2(:,2,:,:,:);
            uV2(:,end,:,:,:) = uV2(:,2,:,:,:);
            uV3(:,:,1,:,:) = uV3(:,:,2,:,:);
            uV3(:,:,end,:,:) = uV3(:,:,end-1,:,:);
            uV4(:,:,:,1,:) = uV4(:,:,:,2,:);
            uV4(:,:,:,end,:) = uV4(:,:,:,end-1,:);
            uA(:,:,:,:,1) = uA(:,:,:,:,2);
            uA(:,:,:,:,end) = uA(:,:,:,:,end-1);
            UV{t} = {uV1,uV2,uV3,uV4,uA};
            
            % ==================== Finite Differences Section ====================>>

            fprintf('Delta T / CFL = %.8f.\n',dt/min(min(min(min(min(CFL))))));
            
            CFL_Plot(t,:) = [dt/min(min(min(min(min(CFL))))),max(max(max(max(max(CFL_v1))))),...
                max(max(max(max(max(CFL_v2))))),max(max(max(max(max(CFL_v3))))),max(max(max(max(max(CFL_v4))))),max(max(max(max(max(CFL_A)))))];
            
            % ==================== Controls Section ====================>>

            Size_Mat = size(u);
            uV1_PF = uV1(:);
            uV2_PF = uV2(:);
            uV3_PF = uV3(:);
            uV4_PF = uV4(:);
            uA_PF = uA(:);
            Controls_Aux = zeros(TL,14);
            Warm_Start = reshape(Controls_T,LV1,LV2,LV3,LV4,LA,14);
            Warm_Start = Warm_Start(:,:,:,:,:,1:12);

            if WarmStart == 0
                Warm_Start(Warm_Start ~= 0) = 0;
            elseif WarmStart ~= 1
                disp('Choose WarmStart 0 or 1.');
                return;
            end

            parfor index = 1:TL

                v3_aux = uV3_PF(index);
                v2_aux = uV2_PF(index);
                
                if t >= Z32
                    if olCV3(v3_aux,t) > 0
                        Aux_V3 = 0;
                    else
                        Aux_V3 = 1;
                    end
                    if artificialVirtualControls
                        Aux_V3 = 0;
                    end
                else
                    Aux_V3 = Phi_V3(t);
                end
                
                if t >= Z21
                    if olCV2(v2_aux,t) > 0
                        Aux_V2 = 0;
                    else
                        Aux_V2 = 1;
                    end
                    if artificialVirtualControls
                        Aux_V2 = 0;
                    end
                else
                    Aux_V2 = Phi_V2(t);
                end

                [ind_V1,ind_V2,ind_V3,ind_V4,ind_A] = ind2sub(Size_Mat,index);
                WS = squeeze(Warm_Start(ind_V1,ind_V2,ind_V3,ind_V4,ind_A,:));
                Aux_uA = uA_PF(index);
                Aux_Lims = limsA(A(ind_A));

                Controls_Aux(index,:) = [Quad_FMC_Lambda_8(Q(V1(ind_V1),V2(ind_V2),V3(ind_V3),t),b(V1(ind_V1),V2(ind_V2),V3(ind_V3),...
                        V4(ind_V4),t),c(t),d(V1(ind_V1),V2(ind_V2),V3(ind_V3),uV1_PF(index),v2_aux,v3_aux,uV4_PF(index),...
                        Aux_uA,t),dd,options,WS,Aux_Lims,InfTermicas);Aux_V2;Aux_V3];

            end

            Controls_T = reshape(Controls_Aux,LV1,LV2,LV3,LV4,LA,14);

            % ==================== Time Step ====================>>

            if t ~= length(time)
                for v1 = 1:length(V1)
                    for v2 = 1:length(V2)
                        for v3 = 1:length(V3)
                            for v4 = 1:length(V4)
                                for a = 1:length(A)

                                    Controls(1:14) = ControlsFMC_Aux(v1,v2,v3,v4,a,:);

                                    HAM(v1,v2,v3,v4,a) = TMax*(Controls*[d(V1(v1),V2(v2),V3(v3),uV1(v1,v2,v3,v4,a),...
                                        uV2(v1,v2,v3,v4,a),uV3(v1,v2,v3,v4,a),uV4(v1,v2,v3,v4,a),uA(v1,v2,v3,v4,a),t);olCV2(uV2(v1,v2,v3,v4,a),t);...
                                        olCV3(uV3(v1,v2,v3,v4,a),t)]+H0(uV1(v1,v2,v3,v4,a),uV1(v1,v2,v3,v4,a),uV3(v1,v2,v3,v4,a),uV4(v1,v2,v3,v4,a),t));

                                    u(v1,v2,v3,v4,a) = u(v1,v2,v3,v4,a) + dt*HAM(v1,v2,v3,v4,a);

                                end
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
                    surf(vplot3_1,vplot4_1,squeeze(u(index,1,:,:,1)/NormCost)');
                    title(['Optimal Cost (A = 0) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',3);
                clf(3); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(u(index,1,:,:,end)/NormCost)');
                    title(['Optimal Cost (A = 1) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',4);
                clf(4); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(u(index,1,:,:,1+floor(length(A)/2))/NormCost)');
                    title(['Optimal Cost (A = 0.5) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',5);
                clf(5); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(HAM(index,1,:,:,1)/NormCost)');
                    title(['Hamiltonian (A = 0) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',6);
                clf(6); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(HAM(index,1,:,:,end)/NormCost)');
                    title(['Hamiltonian (A = 1) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',7);
                clf(7); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(HAM(index,1,:,:,1+floor(length(A)/2))/NormCost)');
                    title(['Hamiltonian (A = 0.5) at time-step = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',8);
                clf(8); hold on;
                for i = 1:floor(length(V3)/2):length(V3)
                    surf(vplot4_A,vplotA_4,squeeze(u(end,1,i,:,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['Optimal Cost ($V1 = \underline{V}_1$) at time = ',num2str(t)],'Interpreter','latex');
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',9);
                clf(9); hold on;
                for index = 1:floor(length(V4)/2):length(V4)
                    surf(vplot3_A,vplotA_3,squeeze(uV3(index,1,:,index,:)/NormCost));
                end
                xlabel('Vol. Palmar'); ylabel('Charge Battery');
                title(['u_{V3} at time = ',num2str(t)]);
                xlim([Vmin3 Vmax3]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',10);
                clf(10); hold on;
                for index = 1:length(V4)
                    surf(vplot3_A,vplotA_3,squeeze(uV1(index,1,:,index,:)/NormCost));
                end
                xlabel('Vol. Palmar'); ylabel('Charge Battery');
                title(['u_{V1} at time = ',num2str(t)]);
                xlim([Vmin3 Vmax3]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',14);
                clf(14); hold on;
                for i = 1:floor(length(V4)/2):length(V4)
                    surf(vplot4_A,vplotA_4,squeeze(u(i,1,i,:,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['Optimal Cost at time = ',num2str(t)]);
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',12);
                clf(12); hold on;
                for index = 1:length(V4)
                    surf(vplot4_A,vplotA_4,squeeze(uA(index,1,:,index,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['u_{A} at time = ',num2str(t)]);
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);
                
                set(0,'CurrentFigure',13);
                clf(13); hold on;
                Vols = linspace(min([Vmin1,Vmin3,Vmin4]),max([Vmax1,Vmax3,Vmax4]),1000);
                plot(Vols,K11(Vols,Vols)+K21(Vols,Vols),Vols,K12(Vols,Vols)+K22(Vols,Vols),...
                    Vols,K13(Vols,t)+K23(Vols,t),Vols,K14(Vols,t)+K24);
                P = plot(Vols,K12(Vols,Vols)+K22(Vols,Vols)+K13(Vols,t)+K23(Vols,t));
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
        IP2 = min(abs(V2-Vini2)) == abs(V2-Vini2);
        IP3 = min(abs(V3-Vini3)) == abs(V3-Vini3);
        IP4 = min(abs(V4-Vini4)) == abs(V4-Vini4);
        IPA = min(abs(A-A0)) == abs(A-A0);
        if sum(IP1) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP1)
                    if IP1(r) == 1 && Aux_Cont == 0
                        IP1(r) = 1;
                        Aux_Cont = 1;
                    elseif IP1(r) == 1 && Aux_Cont == 1
                        IP1(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IP2) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP2)
                    if IP2(r) == 1 && Aux_Cont == 0
                        IP2(r) = 1;
                        Aux_Cont = 1;
                    elseif IP2(r) == 1 && Aux_Cont == 1
                        IP2(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IP3) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP3)
                    if IP3(r) == 1 && Aux_Cont == 0
                        IP3(r) = 1;
                        Aux_Cont = 1;
                    elseif IP3(r) == 1 && Aux_Cont == 1
                        IP3(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IP4) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP4)
                    if IP4(r) == 1 && Aux_Cont == 0
                        IP4(r) = 1;
                        Aux_Cont = 1;
                    elseif IP4(r) == 1 && Aux_Cont == 1
                        IP4(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IPA) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IPA)
                    if IPA(r) == 1 && Aux_Cont == 0
                        IPA(r) = 1;
                        Aux_Cont = 1;
                    elseif IPA(r) == 1 && Aux_Cont == 1
                        IPA(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
        Final_Cost = u(IP1,IP2,IP3,IP4,IPA);

        ToSave = {Final_Cost,UV,CFL_Plot,U};
        save([pwd '/','Simulations/Simulation_',name,'.mat'],'ToSave');
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
            
        load([pwd '/','Simulations/Simulation_',name,'.mat']);
        
        if WhatToDo == 7
            
            % We load all the real production which is in MW. We pass it to
            % kw, and we extend the number of points from 145 to
            % length(time).
            
            histPowerBon = Matrix{1}(:,2)*1e3;
            histPowerBay = Matrix{1}(:,3)*1e3;
            histPowerPal = Matrix{1}(:,4)*1e3;
            histPowerSG = Matrix{1}(:,1)*1e3;
            histPowerMot = Matrix{1}(:,70)*1e3;
            histPowerPTA = Matrix{1}(:,69)*1e3;
            histPowerPTB = Matrix{1}(:,72)*1e3;
            histPowerCTR = Matrix{1}(:,71)*1e3;
            
            realDataTime = linspace(0,1,145); % length(histPowerBon) = 145;
            
            histPowerBon = interp1(realDataTime,histPowerBon,time);
            histPowerBay = interp1(realDataTime,histPowerBay,time);
            histPowerPal = interp1(realDataTime,histPowerPal,time);
            histPowerSG = interp1(realDataTime,histPowerSG,time);
            histPowerMot = interp1(realDataTime,histPowerMot,time);
            histPowerPTA = interp1(realDataTime,histPowerPTA,time);
            histPowerPTB = interp1(realDataTime,histPowerPTB,time);
            histPowerCTR = interp1(realDataTime,histPowerCTR,time);

        end
        
        Final_Cost = ToSave{1};
        if length(Final_Cost) > 1
            Final_Cost(2:end) = [];
        end
        UV = ToSave{2};
        U = ToSave{4};
        CFL_Plot = ToSave{3};
        
        if WhatToDo ~= 5
        
            IP1 = min(abs(V1-Vini1)) == abs(V1-Vini1);
            IP2 = min(abs(V2-Vini2)) == abs(V2-Vini2);
            IP3 = min(abs(V3-Vini3)) == abs(V3-Vini3);
            IP4 = min(abs(V4-Vini4)) == abs(V4-Vini4);
            IPA = min(abs(A-A0)) == abs(A-A0);

            Optimal_Cost = zeros(1,length(time));
            instantTotalCost = zeros(1,length(time));
            futureCostOP = zeros(1,length(time));
            futureRealCostOP = zeros(1,length(time));
            Phi_V2_Aux = zeros(1,length(time)); % To save the past virtual control V2.
            Phi_V3_Aux = zeros(1,length(time)); % To save the past virtual control V3.
            Power_H1_OP = zeros(1,length(time));
            Power_H2_OP = zeros(1,length(time));
            Power_H3_OP = zeros(1,length(time));
            Power_H4_OP = zeros(1,length(time));
            Vol_OP_1 = zeros(1,length(time));
            Vol_OP_2 = zeros(1,length(time));
            Vol_OP_3 = zeros(1,length(time));
            Vol_OP_4 = zeros(1,length(time));
            Vol_OP_Real_1 = zeros(1,length(time));
            Vol_OP_Real_2 = zeros(1,length(time));
            Vol_OP_Real_3 = zeros(1,length(time));
            Vol_OP_Real_4 = zeros(1,length(time));
            if sum(IP1) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP1)
                    if IP1(r) == 1 && Aux_Cont == 0
                        IP1(r) = 1;
                        Aux_Cont = 1;
                    elseif IP1(r) == 1 && Aux_Cont == 1
                        IP1(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IP2) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP2)
                    if IP2(r) == 1 && Aux_Cont == 0
                        IP2(r) = 1;
                        Aux_Cont = 1;
                    elseif IP2(r) == 1 && Aux_Cont == 1
                        IP2(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IP3) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP3)
                    if IP3(r) == 1 && Aux_Cont == 0
                        IP3(r) = 1;
                        Aux_Cont = 1;
                    elseif IP3(r) == 1 && Aux_Cont == 1
                        IP3(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            if sum(IP4) ~= 1
                Aux_Cont = 0;
                for r = 1:length(IP4)
                    if IP4(r) == 1 && Aux_Cont == 0
                        IP4(r) = 1;
                        Aux_Cont = 1;
                    elseif IP4(r) == 1 && Aux_Cont == 1
                        IP4(r) = 0;
                        Aux_Cont = 1;
                    end
                end
            end
            Vol_OP_1(1) = V1(IP1);
            Vol_OP_2(1) = V2(IP2);
            Vol_OP_3(1) = V3(IP3);
            Vol_OP_4(1) = V4(IP4);
            Vol_OP_Real_1(1) = V1(IP1);
            Vol_OP_Real_2(1) = V2(IP2);
            Vol_OP_Real_3(1) = V3(IP3);
            Vol_OP_Real_4(1) = V4(IP4);
            Spot = zeros(1,length(time));
            costControllable = zeros(12,length(time));
            totalInstantCost = zeros(4,length(time));
            NC_Cost = zeros(1,length(time));
            uV1plot = zeros(1,length(time));
            uV2plot = zeros(1,length(time));
            uV3plot = zeros(1,length(time));
            uV4plot = zeros(1,length(time));
            uPlot = zeros(1,length(time));
            accumContCost = zeros(11,length(time));
            accumRealInstantCost = zeros(1,length(time));
            accumCost = zeros(11,length(time));
            instantCostHJB =  zeros(1,length(time));
            realInstantCost = zeros(1,length(time));
            WaterBonete           = zeros(1,length(time)); 
            WaterBaygorria        = zeros(1,length(time)); 
            WaterPalmar           = zeros(1,length(time)); 
            WaterSG               = zeros(1,length(time)); 
            WaterVirtualBaygorria = zeros(1,length(time));
            WaterVirtualPalmar    = zeros(1,length(time));
            dPH1overFlowPlot = zeros(1,length(time));
            dPH2overFlowPlot = zeros(1,length(time));
            spotTimesPH1 = zeros(1,length(time));
            spotTimesPH2 = zeros(1,length(time));
            CV2InstantCost = zeros(1,length(time));
            CV3InstantCost = zeros(1,length(time));
            
            waterBoneteOverTime      = zeros(1,length(time));
            waterBaygorriaOverTime   = zeros(1,length(time));
            waterPalmarOverTime      = zeros(1,length(time));
            waterSGOverTime          = zeros(1,length(time));
            virtualBaygorriaOverTime = zeros(1,length(time));
            virtualPalmarOverTime    = zeros(1,length(time));
            
            NC_Cost_Acum = zeros(1,length(time));
            
            if WhatToDo ~= 7
                Controls_OP = zeros(length(time),14);
            end
            
            CostControls = zeros(12,length(time));
            Battery_Real = zeros(1,length(time));
            Power_Battery_OP = zeros(1,length(time));
            Prec_Cont = zeros(length(time)+1,12);
            Power_F1_OP = zeros(1,length(time));
            Power_F2_OP = zeros(1,length(time));
            Power_F3_OP = zeros(1,length(time));
            Power_F4_OP = zeros(1,length(time));
            Power_F5_OP = zeros(1,length(time));
            Battery_Cost = 0;
            lambda = {};

            uAplot = zeros(1,length(time));
            Battery_OP = A(IPA);
            Battery_Real(1) = A(IPA);
            x0 = zeros(12,1);

            Power_H1 = @(x1,x2,V1,V2) eta1*x1*FlowMax1(V1,V2)*(dH1(V1,V2)-d1(V1,V2,x1,x2)); % In kW.
            Power_H2 = @(x3,x4,V2,V3) eta2*x3*FlowMax2(V2,V3)*(dH2(V2,V3)-d2(V2,V3,x3,x4)); % In kW.
            Power_H3 = @(x5,V3,t)     eta3*x5*FlowMax3(V3,t)*(dH3(V3,t)-d3(V3,x5,0,t)); % In kW.
            Power_H4 = @(x6,V4,t)     eta4*x6*FlowMax4*(dH4(V4,t)-d4(x6,0)); % In kW.
            Power_Battery = @(x7)     x7*PAMax; % In kW.
                        
            dPH1overFlow = @(x1,x2,V1,V2) eta1(dH1(V1,V2)-c1*(2*x1*FlowMax1(V1,V2)+x2*MaxSpill1(V1,V2))); % In kW/(m^3/s).
            dPH2overFlow = @(x3,x4,V2,V3) eta2(dH2(V2,V3)-c2*(2*x3*FlowMax2(V2,V3)+x4*MaxSpill2(V2,V3))); % In kW/(m^3/s).
            
            contBon = @(potBon,spillBon,V1,V2) ((c1*spillBon*FlowMax1(V1,V2)*MaxSpill1(V1,V2) - eta1*FlowMax1(V1,V2)*dH1(V1,V2)) +...
                sqrt((c1*spillBon*FlowMax1(V1,V2)*MaxSpill1(V1,V2) - eta1*FlowMax1(V1,V2)*dH1(V1,V2))^2 - 4*(-c1*eta1*FlowMax1(V1,V2))*(-potBon)))/...
                (2*(-c1*eta1*FlowMax1(V1,V2)));
            contBay = @(potBay,spillBay,V2,V3) ((c2*spillBay*FlowMax2(V2,V3)*MaxSpill2(V2,V3) - eta2*FlowMax2(V2,V3)*dH2(V2,V3)) +...
                sqrt((c2*spillBay*FlowMax2(V2,V3)*MaxSpill2(V2,V3) - eta2*FlowMax2(V2,V3)*dH2(V2,V3))^2 - 4*(-c2*eta2*FlowMax2(V2,V3))*(-potBay)))/...
                (2*(-c2*eta2*FlowMax2(V2,V3)));
            contPal = @(potPal,spillPal,V3,t) ((c3*spillPal*FlowMax3(V3,t)*MaxSpill3(V3,t) - eta3*FlowMax3(V3,t)*dH3(V3,t)) +...
                sqrt((c3*spillPal*FlowMax3(V3,t)*MaxSpill3(V3,t) - eta3*FlowMax3(V3,t)*dH3(V3,t))^2 - 4*(-c3*eta3*FlowMax3(V3,t))*(-potPal)))/...
                (2*(-c3*eta3*FlowMax3(V3,t)));
            contSG = @(potSG,spillSG,V4,t) ((c4*spillSG*FlowMax4*MaxSpill4 - eta4*FlowMax4*dH4(V4,t)) +...
                sqrt((c4*spillSG*FlowMax4*MaxSpill4 - eta4*FlowMax4*dH4(V4,t))^2 - 4*(-c4*eta4*FlowMax4)*(-potSG)))/...
                (2*(-c4*eta4*FlowMax4));
            contMot = @(potMot) potMot/MaxTermicas(1);
            contPTA = @(potPTA) potPTA/MaxTermicas(2);
            contPTB = @(potPTB) potPTB/MaxTermicas(3);
            contCTR = @(potCTR) potCTR/MaxTermicas(4);
            contFail = @(potFail) potFail/MaxTermicas(5);            
            
            for t = 1:length(time)

                fprintf('Completed: %.2f.\n',100*(t/length(time)));

                if t == 1
                    uV1 = UV{t}{1}(IP1,IP2,IP3,IP4,IPA);
                    uV2 = UV{t}{2}(IP1,IP2,IP3,IP4,IPA);
                    uV3 = UV{t}{3}(IP1,IP2,IP3,IP4,IPA);
                    uV4 = UV{t}{4}(IP1,IP2,IP3,IP4,IPA);
                    uA = UV{t}{5}(IP1,IP2,IP3,IP4,IPA);
                    u = U{t}(IP1,IP2,IP3,IP4,IPA);
                else
                    realxV1 = (Vol_OP_Real_1(t)-V1(IP1))/(V1(IP1+1)-V1(IP1));
                    realxV2 = (Vol_OP_Real_2(t)-V2(IP2))/(V2(IP2+1)-V2(IP2));
                    realxV3 = (Vol_OP_Real_3(t)-V3(IP3))/(V3(IP3+1)-V3(IP3));
                    realxV4 = (Vol_OP_Real_4(t)-V4(IP4))/(V4(IP4+1)-V4(IP4));
                    realxA = (Battery_Real(t)-A(IPA))/(A(IPA+1)-A(IPA));
                    real = [realxV1,realxV2,realxV3,realxV4,realxA];

                    valsV1 = [UV{t}{1}(IP1,IP2,IP3,IP4,IPA); UV{t}{1}(IP1+1,IP2,IP3,IP4,IPA); UV{t}{1}(IP1,IP2+1,IP3,IP4,IPA); UV{t}{1}(IP1+1,IP2+1,IP3,IP4,IPA);...
                        UV{t}{1}(IP1,IP2,IP3+1,IP4,IPA); UV{t}{1}(IP1+1,IP2,IP3+1,IP4,IPA); UV{t}{1}(IP1,IP2+1,IP3+1,IP4,IPA); UV{t}{1}(IP1+1,IP2+1,IP3+1,IP4,IPA);...
                        UV{t}{1}(IP1,IP2,IP3,IP4+1,IPA); UV{t}{1}(IP1+1,IP2,IP3,IP4+1,IPA); UV{t}{1}(IP1,IP2+1,IP3,IP4+1,IPA); UV{t}{1}(IP1+1,IP2+1,IP3,IP4+1,IPA);...
                        UV{t}{1}(IP1,IP2,IP3+1,IP4+1,IPA); UV{t}{1}(IP1+1,IP2,IP3+1,IP4+1,IPA); UV{t}{1}(IP1,IP2+1,IP3+1,IP4+1,IPA); UV{t}{1}(IP1+1,IP2+1,IP3+1,IP4+1,IPA);...
                        UV{t}{1}(IP1,IP2,IP3,IP4,IPA+1); UV{t}{1}(IP1+1,IP2,IP3,IP4,IPA+1); UV{t}{1}(IP1,IP2+1,IP3,IP4,IPA+1); UV{t}{1}(IP1+1,IP2+1,IP3,IP4,IPA+1);...
                        UV{t}{1}(IP1,IP2,IP3+1,IP4,IPA+1); UV{t}{1}(IP1+1,IP2,IP3+1,IP4,IPA+1); UV{t}{1}(IP1,IP2+1,IP3+1,IP4,IPA+1); UV{t}{1}(IP1+1,IP2+1,IP3+1,IP4,IPA+1);...
                        UV{t}{1}(IP1,IP2,IP3,IP4+1,IPA+1); UV{t}{1}(IP1+1,IP2,IP3,IP4+1,IPA+1); UV{t}{1}(IP1,IP2+1,IP3,IP4+1,IPA+1); UV{t}{1}(IP1+1,IP2+1,IP3,IP4+1,IPA+1);...
                        UV{t}{1}(IP1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{1}(IP1+1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{1}(IP1,IP2+1,IP3+1,IP4+1,IPA+1); UV{t}{1}(IP1+1,IP2+1,IP3+1,IP4+1,IPA+1)];
                    uV1 = Int5D_2(real,valsV1);
                    
                    valsV2 = [UV{t}{2}(IP1,IP2,IP3,IP4,IPA); UV{t}{2}(IP1+1,IP2,IP3,IP4,IPA); UV{t}{2}(IP1,IP2+1,IP3,IP4,IPA); UV{t}{2}(IP1+1,IP2+1,IP3,IP4,IPA);...
                        UV{t}{2}(IP1,IP2,IP3+1,IP4,IPA); UV{t}{2}(IP1+1,IP2,IP3+1,IP4,IPA); UV{t}{2}(IP1,IP2+1,IP3+1,IP4,IPA); UV{t}{2}(IP1+1,IP2+1,IP3+1,IP4,IPA);...
                        UV{t}{2}(IP1,IP2,IP3,IP4+1,IPA); UV{t}{2}(IP1+1,IP2,IP3,IP4+1,IPA); UV{t}{2}(IP1,IP2+1,IP3,IP4+1,IPA); UV{t}{2}(IP1+1,IP2+1,IP3,IP4+1,IPA);...
                        UV{t}{2}(IP1,IP2,IP3+1,IP4+1,IPA); UV{t}{2}(IP1+1,IP2,IP3+1,IP4+1,IPA); UV{t}{2}(IP1,IP2+1,IP3+1,IP4+1,IPA); UV{t}{2}(IP1+1,IP2+1,IP3+1,IP4+1,IPA);...
                        UV{t}{2}(IP1,IP2,IP3,IP4,IPA+1); UV{t}{2}(IP1+1,IP2,IP3,IP4,IPA+1); UV{t}{2}(IP1,IP2+1,IP3,IP4,IPA+1); UV{t}{2}(IP1+1,IP2+1,IP3,IP4,IPA+1);...
                        UV{t}{2}(IP1,IP2,IP3+1,IP4,IPA+1); UV{t}{2}(IP1+1,IP2,IP3+1,IP4,IPA+1); UV{t}{2}(IP1,IP2+1,IP3+1,IP4,IPA+1); UV{t}{2}(IP1+1,IP2+1,IP3+1,IP4,IPA+1);...
                        UV{t}{2}(IP1,IP2,IP3,IP4+1,IPA+1); UV{t}{2}(IP1+1,IP2,IP3,IP4+1,IPA+1); UV{t}{2}(IP1,IP2+1,IP3,IP4+1,IPA+1); UV{t}{2}(IP1+1,IP2+1,IP3,IP4+1,IPA+1);...
                        UV{t}{2}(IP1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{2}(IP1+1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{2}(IP1,IP2+1,IP3+1,IP4+1,IPA+1); UV{t}{2}(IP1+1,IP2+1,IP3+1,IP4+1,IPA+1)];
                    uV2 = Int5D_2(real,valsV2);
                    
                    valsV3 = [UV{t}{3}(IP1,IP2,IP3,IP4,IPA); UV{t}{3}(IP1+1,IP2,IP3,IP4,IPA); UV{t}{3}(IP1,IP2+1,IP3,IP4,IPA); UV{t}{3}(IP1+1,IP2+1,IP3,IP4,IPA);...
                        UV{t}{3}(IP1,IP2,IP3+1,IP4,IPA); UV{t}{3}(IP1+1,IP2,IP3+1,IP4,IPA); UV{t}{3}(IP1,IP2+1,IP3+1,IP4,IPA); UV{t}{3}(IP1+1,IP2+1,IP3+1,IP4,IPA);...
                        UV{t}{3}(IP1,IP2,IP3,IP4+1,IPA); UV{t}{3}(IP1+1,IP2,IP3,IP4+1,IPA); UV{t}{3}(IP1,IP2+1,IP3,IP4+1,IPA); UV{t}{3}(IP1+1,IP2+1,IP3,IP4+1,IPA);...
                        UV{t}{3}(IP1,IP2,IP3+1,IP4+1,IPA); UV{t}{3}(IP1+1,IP2,IP3+1,IP4+1,IPA); UV{t}{3}(IP1,IP2+1,IP3+1,IP4+1,IPA); UV{t}{3}(IP1+1,IP2+1,IP3+1,IP4+1,IPA);...
                        UV{t}{3}(IP1,IP2,IP3,IP4,IPA+1); UV{t}{3}(IP1+1,IP2,IP3,IP4,IPA+1); UV{t}{3}(IP1,IP2+1,IP3,IP4,IPA+1); UV{t}{3}(IP1+1,IP2+1,IP3,IP4,IPA+1);...
                        UV{t}{3}(IP1,IP2,IP3+1,IP4,IPA+1); UV{t}{3}(IP1+1,IP2,IP3+1,IP4,IPA+1); UV{t}{3}(IP1,IP2+1,IP3+1,IP4,IPA+1); UV{t}{3}(IP1+1,IP2+1,IP3+1,IP4,IPA+1);...
                        UV{t}{3}(IP1,IP2,IP3,IP4+1,IPA+1); UV{t}{3}(IP1+1,IP2,IP3,IP4+1,IPA+1); UV{t}{3}(IP1,IP2+1,IP3,IP4+1,IPA+1); UV{t}{3}(IP1+1,IP2+1,IP3,IP4+1,IPA+1);...
                        UV{t}{3}(IP1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{3}(IP1+1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{3}(IP1,IP2+1,IP3+1,IP4+1,IPA+1); UV{t}{3}(IP1+1,IP2+1,IP3+1,IP4+1,IPA+1)];
                    uV3 = Int5D_2(real,valsV3);
                    
                    valsV4 = [UV{t}{4}(IP1,IP2,IP3,IP4,IPA); UV{t}{4}(IP1+1,IP2,IP3,IP4,IPA); UV{t}{4}(IP1,IP2+1,IP3,IP4,IPA); UV{t}{4}(IP1+1,IP2+1,IP3,IP4,IPA);...
                        UV{t}{4}(IP1,IP2,IP3+1,IP4,IPA); UV{t}{4}(IP1+1,IP2,IP3+1,IP4,IPA); UV{t}{4}(IP1,IP2+1,IP3+1,IP4,IPA); UV{t}{4}(IP1+1,IP2+1,IP3+1,IP4,IPA);...
                        UV{t}{4}(IP1,IP2,IP3,IP4+1,IPA); UV{t}{4}(IP1+1,IP2,IP3,IP4+1,IPA); UV{t}{4}(IP1,IP2+1,IP3,IP4+1,IPA); UV{t}{4}(IP1+1,IP2+1,IP3,IP4+1,IPA);...
                        UV{t}{4}(IP1,IP2,IP3+1,IP4+1,IPA); UV{t}{4}(IP1+1,IP2,IP3+1,IP4+1,IPA); UV{t}{4}(IP1,IP2+1,IP3+1,IP4+1,IPA); UV{t}{4}(IP1+1,IP2+1,IP3+1,IP4+1,IPA);...
                        UV{t}{4}(IP1,IP2,IP3,IP4,IPA+1); UV{t}{4}(IP1+1,IP2,IP3,IP4,IPA+1); UV{t}{4}(IP1,IP2+1,IP3,IP4,IPA+1); UV{t}{4}(IP1+1,IP2+1,IP3,IP4,IPA+1);...
                        UV{t}{4}(IP1,IP2,IP3+1,IP4,IPA+1); UV{t}{4}(IP1+1,IP2,IP3+1,IP4,IPA+1); UV{t}{4}(IP1,IP2+1,IP3+1,IP4,IPA+1); UV{t}{4}(IP1+1,IP2+1,IP3+1,IP4,IPA+1);...
                        UV{t}{4}(IP1,IP2,IP3,IP4+1,IPA+1); UV{t}{4}(IP1+1,IP2,IP3,IP4+1,IPA+1); UV{t}{4}(IP1,IP2+1,IP3,IP4+1,IPA+1); UV{t}{4}(IP1+1,IP2+1,IP3,IP4+1,IPA+1);...
                        UV{t}{4}(IP1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{4}(IP1+1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{4}(IP1,IP2+1,IP3+1,IP4+1,IPA+1); UV{t}{4}(IP1+1,IP2+1,IP3+1,IP4+1,IPA+1)];
                    uV4 = Int5D_2(real,valsV4);
                    
                    valsA = [UV{5}{4}(IP1,IP2,IP3,IP4,IPA); UV{t}{5}(IP1+1,IP2,IP3,IP4,IPA); UV{t}{5}(IP1,IP2+1,IP3,IP4,IPA); UV{t}{5}(IP1+1,IP2+1,IP3,IP4,IPA);...
                        UV{t}{5}(IP1,IP2,IP3+1,IP4,IPA); UV{t}{5}(IP1+1,IP2,IP3+1,IP4,IPA); UV{t}{5}(IP1,IP2+1,IP3+1,IP4,IPA); UV{t}{5}(IP1+1,IP2+1,IP3+1,IP4,IPA);...
                        UV{t}{5}(IP1,IP2,IP3,IP4+1,IPA); UV{t}{5}(IP1+1,IP2,IP3,IP4+1,IPA); UV{t}{5}(IP1,IP2+1,IP3,IP4+1,IPA); UV{t}{5}(IP1+1,IP2+1,IP3,IP4+1,IPA);...
                        UV{t}{5}(IP1,IP2,IP3+1,IP4+1,IPA); UV{t}{5}(IP1+1,IP2,IP3+1,IP4+1,IPA); UV{t}{5}(IP1,IP2+1,IP3+1,IP4+1,IPA); UV{t}{5}(IP1+1,IP2+1,IP3+1,IP4+1,IPA);...
                        UV{t}{5}(IP1,IP2,IP3,IP4,IPA+1); UV{t}{5}(IP1+1,IP2,IP3,IP4,IPA+1); UV{t}{5}(IP1,IP2+1,IP3,IP4,IPA+1); UV{t}{5}(IP1+1,IP2+1,IP3,IP4,IPA+1);...
                        UV{t}{5}(IP1,IP2,IP3+1,IP4,IPA+1); UV{t}{5}(IP1+1,IP2,IP3+1,IP4,IPA+1); UV{t}{5}(IP1,IP2+1,IP3+1,IP4,IPA+1); UV{t}{5}(IP1+1,IP2+1,IP3+1,IP4,IPA+1);...
                        UV{t}{5}(IP1,IP2,IP3,IP4+1,IPA+1); UV{t}{5}(IP1+1,IP2,IP3,IP4+1,IPA+1); UV{t}{5}(IP1,IP2+1,IP3,IP4+1,IPA+1); UV{t}{5}(IP1+1,IP2+1,IP3,IP4+1,IPA+1);...
                        UV{t}{5}(IP1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{5}(IP1+1,IP2,IP3+1,IP4+1,IPA+1); UV{t}{5}(IP1,IP2+1,IP3+1,IP4+1,IPA+1); UV{t}{5}(IP1+1,IP2+1,IP3+1,IP4+1,IPA+1)];
                    uA = Int5D_2(real,valsA);
                    
                    valsU = [U{t}(IP1,IP2,IP3,IP4,IPA); U{t}(IP1+1,IP2,IP3,IP4,IPA); U{t}(IP1,IP2+1,IP3,IP4,IPA); U{t}(IP1+1,IP2+1,IP3,IP4,IPA);...
                        U{t}(IP1,IP2,IP3+1,IP4,IPA); U{t}(IP1+1,IP2,IP3+1,IP4,IPA); U{t}(IP1,IP2+1,IP3+1,IP4,IPA); U{t}(IP1+1,IP2+1,IP3+1,IP4,IPA);...
                        U{t}(IP1,IP2,IP3,IP4+1,IPA); U{t}(IP1+1,IP2,IP3,IP4+1,IPA); U{t}(IP1,IP2+1,IP3,IP4+1,IPA); U{t}(IP1+1,IP2+1,IP3,IP4+1,IPA);...
                        U{t}(IP1,IP2,IP3+1,IP4+1,IPA); U{t}(IP1+1,IP2,IP3+1,IP4+1,IPA); U{t}(IP1,IP2+1,IP3+1,IP4+1,IPA); U{t}(IP1+1,IP2+1,IP3+1,IP4+1,IPA);...
                        U{t}(IP1,IP2,IP3,IP4,IPA+1); U{t}(IP1+1,IP2,IP3,IP4,IPA+1); U{t}(IP1,IP2+1,IP3,IP4,IPA+1); U{t}(IP1+1,IP2+1,IP3,IP4,IPA+1);...
                        U{t}(IP1,IP2,IP3+1,IP4,IPA+1); U{t}(IP1+1,IP2,IP3+1,IP4,IPA+1); U{t}(IP1,IP2+1,IP3+1,IP4,IPA+1); U{t}(IP1+1,IP2+1,IP3+1,IP4,IPA+1);...
                        U{t}(IP1,IP2,IP3,IP4+1,IPA+1); U{t}(IP1+1,IP2,IP3,IP4+1,IPA+1); U{t}(IP1,IP2+1,IP3,IP4+1,IPA+1); U{t}(IP1+1,IP2+1,IP3,IP4+1,IPA+1);...
                        U{t}(IP1,IP2,IP3+1,IP4+1,IPA+1); U{t}(IP1+1,IP2,IP3+1,IP4+1,IPA+1); U{t}(IP1,IP2+1,IP3+1,IP4+1,IPA+1); U{t}(IP1+1,IP2+1,IP3+1,IP4+1,IPA+1)];
                    u = Int5D_2(real,valsU);
                    
                    if postWithZeroPartials
                        uV1 = 0; uV2 = 0; uV3 = 0; uV4 = 0; uA = 0;
                    end
                    
                end

                ContLims = limsA(Battery_Real(t));
                uV1plot(t) = uV1;
                uV2plot(t) = uV2;
                uV3plot(t) = uV3;
                uV4plot(t) = uV4;
                uAplot(t) = uA;
                uPlot(t) = u;

                % ==================== Controls Section ====================>>

                if t >= Z21
                    if UseLamb21 == 1
                        Real2 = Phi_V2_Aux(t-Z21+1);
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
                    % We do not want to penalize the initial controls.
                    dd_OP = dd;
                    d_OP = @(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t);
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t);
                else
                    dd_OP = dd + diag(rho/dt*ones(1,12)) - diag(rho/dt*ones(1,12).*[0,1,0,0,0,0,0,0,0,0,0,0])*9900/10000;
                    d_OP = @(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t) + ...
                    -2*rho/dt*Conts' + 2*rho/dt*Conts(2)*[0,1,0,0,0,0,0,0,0,0,0,0]';
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t) + rho*sum(Conts.^2)/dt - rho*sum(Conts(2)^2)/dt*[0,1,0,0,0,0,0,0,0,0,0,0]';
                    Added_Cost = @(New_Conts,Conts) sum(New_Conts.^2)*rho/dt - 2*rho/dt*sum(New_Conts.*Conts) + rho/dt*sum(Conts.^2); % <<< I have to correct that, the spillage of Bonete now is L2.
                end

                if WhatToDo ~= 7

                    [Controls_OP(t,1:12)] = Quad_FMC_Lambda_8(Q(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),t),...
                        b(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),Vol_OP_Real_4(t),t),c(t),...
                        d_OP(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),uV1,uV2,uV3,uV4,uA,t,Prec_Cont(t,:)),dd_OP,...
                        options,x0,ContLims,InfTermicas);

                    dd_L = dd;
                    d_L = @(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t) d(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t);
                    % To compute the spot price, I remove the penalization
                    % in the controls, To do this, we redefine dd_L and
                    % d_L, and we use the non-penalized values.
                    
                    [~,lambda{t}] = Quad_FMC_Lambda_9(Q(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),t),...
                        b(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),Vol_OP_Real_4(t),t),c(t),...
                        d_L(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),uV1,uV2,uV3,uV4,uA,t),dd_L,...
                        options,Controls_OP(t,1:12)',ContLims,InfTermicas);

                elseif WhatToDo == 7
                    
                    Controls_OP(t,1) = contBon(histPowerBon(t),0,Vol_OP_Real_1(t),Vol_OP_Real_2(t));
                    Controls_OP(t,2) = 0; % I am assuming that historically, the spillage was zero.
                    Controls_OP(t,3) = contBay(histPowerBay(t),0,Vol_OP_Real_2(t),Vol_OP_Real_3(t));
                    Controls_OP(t,4) = 0; % I am assuming that historically, the spillage was zero.
                    Controls_OP(t,5) = contPal(histPowerPal(t),0,Vol_OP_Real_3(t),t);
                    Controls_OP(t,6) = contPal(histPowerSG(t),0,Vol_OP_Real_4(t),t);
                    Controls_OP(t,7) = 0; % Uruguay does not have a battery.
                    Controls_OP(t,8) = contMot(histPowerMot(t));
                    Controls_OP(t,9) = contPTA(histPowerPTA(t));
                    Controls_OP(t,10) = contPTB(histPowerPTB(t));
                    Controls_OP(t,11) = contCTR(histPowerCTR(t));
                    Controls_OP(t,12) = 0; % I am assuming that Uruguay does not fail.
                    
                end

                Phi_V2_Aux(t) = (Controls_OP(t,1)*FlowMax1(Vol_OP_Real_1(t),Vol_OP_Real_2(t)) + Controls_OP(t,2)*MaxSpill1(Vol_OP_Real_1(t),Vol_OP_Real_2(t))) / MaxV2;
                Phi_V3_Aux(t) = (Controls_OP(t,3)*FlowMax2(Vol_OP_Real_2(t),Vol_OP_Real_3(t)) + Controls_OP(t,4)*MaxSpill2(Vol_OP_Real_2(t),Vol_OP_Real_3(t))) / MaxV3;
                x0(1:12) = Controls_OP(t,1:12);
                Prec_Cont(t+1,1:12) = Controls_OP(t,1:12);
                
                %%%%% This is related to the virtual flows (DOWN).
                if UseLamb32 == 1
                    if t >= Z32
                        Controls_OP(t,14) = Phi_V3_Aux(t-Z32+1);
                        if Fixed == 0
                            if olCV3(uV3,t) > 0
                                Controls_OP(t,14) = 0;
                            else
                                Controls_OP(t,14) = 1;
                            end
                        end
                        if artificialVirtualControls
                            Controls_OP(t,14) = 0;
                        end
                    else
                        Controls_OP(t,14) = Phi_V3(t);
                    end
                else
                    Controls_OP(t,14) = 0;
                end
                
                if UseLamb21 == 1
                    if t >= Z21
                        Controls_OP(t,13) = Phi_V2_Aux(t-Z21+1);
                        if Fixed == 0
                            if olCV2(uV2,t) > 0
                                Controls_OP(t,13) = 0;
                            else
                                Controls_OP(t,13) = 1;
                            end
                        end
                        if artificialVirtualControls
                            Controls_OP(t,13) = 0;
                        end
                    else
                        Controls_OP(t,13) = Phi_V2(t);
                    end
                else
                    Controls_OP(t,13) = 0;
                end
                CV2InstantCost(t) = olCV2(uV2,t);
                CV3InstantCost(t) = olCV3(uV3,t);
                %%%%% This is related to the virtual flows (UP).
                
                % This is the cost vector x^T.d (I did not included x^T.dd.x).
                contInstantCost = dt*TMax*Controls_OP(t,:).*[d(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),uV1,uV2,uV3,uV4,uA,t);olCV2(uV2,t);olCV3(uV3,t)]';
                % This is the cost associated with the non-controllable
                % elements in the Hamiltonian.
                % Notice that: H = min x^T.d + x^T.dd.x + H0.
                NC_Cost(t) = dt*TMax*H0(uV1,uV2,uV3,uV4,t);
                
                % This is the real cost. Here we do not include the partial
                % derivatives, and it is associated with the real optimization
                % problem.
                realInstantCost(t) = dt*TMax*(sum(Controls_OP(t,:).*...
                    [FlowMax1(Vol_OP_Real_1(t),Vol_OP_Real_2(t))*CH1,...
                    MaxSpill1(Vol_OP_Real_1(t),Vol_OP_Real_2(t))*CH1,...
                    FlowMax2(Vol_OP_Real_2(t),Vol_OP_Real_3(t))*CH2,...
                    MaxSpill2(Vol_OP_Real_2(t),Vol_OP_Real_3(t))*CH2,...
                    FlowMax3(Vol_OP_Real_3(t),t)*CH3,...
                    FlowMax4*CH4,...
                    0,... $ Battery.
                    CostTermicasNorm.*MaxTermicas,...
                    MaxV2*(-CH2),...
                    MaxV3*(-CH3)]) - ...
                    (IT1(t)*CH1+IT2(t)*CH2+IT3(t)*CH3+IT4(t)*CH4));
                
                % This is the evaluation of the dual function.
                if WhatToDo == 6 && t >= Z21
                    sft21 = t - Z21 + 1;
                    realInstantCost(t) = realInstantCost(t) + ...
                        dt*TMax*Lambda21(t)*...
                        (Controls_OP(t,13)*MaxV2-Controls_OP(sft21,1)*FlowMax1(Vol_OP_Real_1(sft21),Vol_OP_Real_2(sft21))-Controls_OP(sft21,2)*MaxSpill1(Vol_OP_Real_1(sft21),Vol_OP_Real_2(sft21)));
                    if WhatToDo == 6 && t >= Z32
                    sft32 = t - Z32 + 1;
                    realInstantCost(t) = realInstantCost(t) + ...
                        dt*TMax*Lambda32(t)*...
                        (Controls_OP(t,14)*MaxV3-Controls_OP(sft32,3)*FlowMax2(Vol_OP_Real_2(sft32),Vol_OP_Real_3(sft32))-Controls_OP(sft32,4)*MaxSpill2(Vol_OP_Real_2(sft32),Vol_OP_Real_3(sft32)));
                    end
                end
                
                % These are the controllable costs at time t.
                costControllable(1,t) = contInstantCost(1)+contInstantCost(2); % Controllable cost Bonete.
                costControllable(2,t) = contInstantCost(3)+contInstantCost(4); % Controllable cost Baygorria.
                costControllable(3,t) = contInstantCost(5); % Controllable cost Palmar.
                costControllable(4,t) = contInstantCost(6); % Controllable cost SG.
                costControllable(5,t) = contInstantCost(8); % Controllable cost Motores Batlle.
                costControllable(6,t) = contInstantCost(9); % Controllable cost PTA.
                costControllable(7,t) = contInstantCost(10); % Controllable cost PTB.
                costControllable(8,t) = contInstantCost(11); % Controllable cost CTR.
                costControllable(9,t) = contInstantCost(7); % Controllable cost Battery.
                costControllable(10,t) = contInstantCost(12); % Controllable cost Failure.
                costControllable(11,t) = contInstantCost(13); % Controllable cost Virtual 2.
                costControllable(12,t) = contInstantCost(14); % Controllable cost Virtual 3.
                % The above vector is:
                % [Bonete,Baygorria,Palmar,SG,MB,PTA,PTB,CTR,Battery,Failure].
                
                % These are the total costs at time t.
                totalInstantCost(1,t) = costControllable(1,t) + TMax*dt*H0_1(uV1,t); % Total instant cost Bonete.
                totalInstantCost(2,t) = costControllable(2,t) + TMax*dt*H0_2(uV2,t); % Total instant cost Baygorria.
                totalInstantCost(3,t) = costControllable(3,t) + TMax*dt*H0_3(uV3,t); % Total instant cost Palmar.
                totalInstantCost(4,t) = costControllable(4,t) + TMax*dt*H0_4(uV4,t); % Total instant cost SG.
                                
                % These are the accumulated controllable costs at time t.
%                 if t == 1 ccc
%                     accumRealInstantCost(t) = 0; %realInstantCost(t); ccc
%                 end ccc
                
                if t < length(time)
                    accumRealInstantCost(t+1) = accumRealInstantCost(t) + realInstantCost(t);
                end
                
                if t > 1
%                     accumRealInstantCost(t) = accumRealInstantCost(t-1) + realInstantCost(t); ccc
                    
                    accumContCost(1,t) = accumContCost(1,t-1) + costControllable(1,t); % Accumulated controllable cost Bonete.
                    accumContCost(2,t) = accumContCost(2,t-1) + costControllable(2,t); % Accumulated controllable cost Baygorria.
                    accumContCost(3,t) = accumContCost(3,t-1) + costControllable(3,t); % Accumulated controllable cost Palmar.
                    accumContCost(4,t) = accumContCost(4,t-1) + costControllable(4,t); % Accumulated controllable cost SG.
                    accumContCost(5,t) = accumContCost(5,t-1) + costControllable(5,t); % Accumulated controllable cost Motores Batlle.
                    accumContCost(6,t) = accumContCost(6,t-1) + costControllable(6,t); % Accumulated controllable cost PTA.
                    accumContCost(7,t) = accumContCost(7,t-1) + costControllable(7,t); % Accumulated controllable cost PTB.
                    accumContCost(8,t) = accumContCost(8,t-1) + costControllable(8,t); % Accumulated controllable cost CTR.
                    accumContCost(9,t) = accumContCost(9,t-1) + costControllable(10,t); % Accumulated controllable cost Failure.
                    accumContCost(10,t) = accumContCost(10,t-1) + costControllable(11,t); % Accumulated controllable cost Virtual 2.
                    accumContCost(11,t) = accumContCost(11,t-1) + costControllable(12,t); % Accumulated controllable cost Virtual 3.
                    
                    NC_Cost_Acum(t) = NC_Cost_Acum(t-1) + NC_Cost(t);
                    
                    accumCost(1,t) = accumCost(1,t-1) + totalInstantCost(1,t); % Accumulated total cost Bonete.
                    accumCost(2,t) = accumCost(2,t-1) + totalInstantCost(2,t); % Accumulated total cost Baygorria.
                    accumCost(3,t) = accumCost(3,t-1) + totalInstantCost(3,t); % Accumulated total cost Palmar.
                    accumCost(4,t) = accumCost(4,t-1) + totalInstantCost(4,t); % Accumulated total cost SG.
                    accumCost(5,t) = accumCost(5,t-1) + costControllable(5,t); % Accumulated total cost MB.
                    accumCost(6,t) = accumCost(6,t-1) + costControllable(6,t); % Accumulated total cost PTA.
                    accumCost(7,t) = accumCost(7,t-1) + costControllable(7,t); % Accumulated total cost PTB.
                    accumCost(8,t) = accumCost(8,t-1) + costControllable(8,t); % Accumulated total cost CTR.
                    accumCost(9,t) = accumCost(9,t-1) + costControllable(10,t); % Accumulated total cost Failure.
                    accumCost(10,t) = accumCost(10,t-1) + costControllable(11,t); % Accumulated total cost Virtual 2.
                    accumCost(11,t) = accumCost(11,t-1) + costControllable(12,t); % Accumulated total cost Virtual 3.
                end
                
                CostControls(:,t) = d(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),uV1,uV2,uV3,uV4,uA,t)./...
                   [Power_H1(1,Controls_OP(t,2),Vol_OP_Real_1(t),Vol_OP_Real_2(t)),MaxSpill1(Vol_OP_Real_1(t),Vol_OP_Real_2(t))*dt*TMax/m3Tohm3,...
                   Power_H2(1,Controls_OP(t,4),Vol_OP_Real_2(t),Vol_OP_Real_3(t)),MaxSpill2(Vol_OP_Real_2(t),Vol_OP_Real_3(t))*dt*TMax/m3Tohm3,...
                   Power_H3(1,Vol_OP_Real_3(t),t),Power_H4(1,Vol_OP_Real_4(t),t),PAMax,MaxTermicas]';
               % This cost is in USD/kW for all but the spillages, and is the cost when the devices are working
               % at their maximum power. For the spillages is USD/Control.

                % ==================== Dynamics and Powers Section ====================>>

                if t ~= length(time)

                    Vol_OP_Real_1(t+1) = Vol_OP_Real_1(t) + (dt*TMax/VolMax1)*(IT1(t)-...
                        Controls_OP(t,1)*FlowMax1(Vol_OP_1(t),Vol_OP_2(t))-Controls_OP(t,2)*MaxSpill1(Vol_OP_1(t),Vol_OP_2(t)));
                    Vol_OP_Real_2(t+1) = Vol_OP_Real_2(t) + (dt*TMax/VolMax1)*(IT2(t)-...
                        Controls_OP(t,3)*FlowMax2(Vol_OP_2(t),Vol_OP_3(t))-Controls_OP(t,4)*MaxSpill2(Vol_OP_2(t),Vol_OP_3(t))+Controls_OP(t,13)*MaxV2);
                    Vol_OP_Real_3(t+1) = Vol_OP_Real_3(t) + (dt*TMax/VolMax3)*(IT3(t)-...
                        Controls_OP(t,5)*FlowMax3(Vol_OP_3(t),t)+Controls_OP(t,14)*MaxV3);
                    Vol_OP_Real_4(t+1) = Vol_OP_Real_4(t) + (dt*TMax/VolMax4)*(IT4(t)-...
                        Controls_OP(t,6)*FlowMax4);
                    Battery_Real(t+1) = Battery_Real(t) - dt*TMax*K_A*Controls_OP(t,7);

                    AuxV1 = V1(min(abs(V1-Vol_OP_Real_1(t+1))) == abs(V1-Vol_OP_Real_1(t+1)));
                    Vol_OP_1(t+1) = AuxV1(1);
                    AuxV2 = V2(min(abs(V2-Vol_OP_Real_2(t+1))) == abs(V2-Vol_OP_Real_2(t+1)));
                    Vol_OP_2(t+1) = AuxV2(1);
                    AuxV3 = V3(min(abs(V3-Vol_OP_Real_3(t+1))) == abs(V3-Vol_OP_Real_3(t+1)));
                    Vol_OP_3(t+1) = AuxV3(1);
                    AuxV4 = V4(min(abs(V4-Vol_OP_Real_4(t+1))) == abs(V4-Vol_OP_Real_4(t+1)));
                    Vol_OP_4(t+1) = AuxV4(1);
                    AuxA = A(min(abs(A-Battery_Real(t+1))) == abs(A-Battery_Real(t+1)));
                    Battery_OP(t+1) = AuxA(1);

                    IP1 = find(V1 == Vol_OP_1(t+1));
                    IP2 = find(V2 == Vol_OP_2(t+1));
                    IP3 = find(V3 == Vol_OP_3(t+1));
                    IP4 = find(V4 == Vol_OP_4(t+1));
                    IPA = find(A == Battery_OP(t+1));

                    %%%%% This fixes the volume of the dams between what we
                    % define as maximum and minimum during the simulation
                    % of the HJB equation (DOWN).
                    if Vol_OP_1(t+1) > Vol_OP_Real_1(t+1)
                        IP1 = IP1 - 1;
                    end
                    if Vol_OP_2(t+1) > Vol_OP_Real_2(t+1)
                        IP2 = IP2 - 1;
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
                    
                    if Vol_OP_Real_2(t+1) >= Vmax2
                        Vol_OP_Real_2(t+1) = Vmax2;
                        IP2 = length(V2)-1;
                    elseif Vol_OP_Real_2(t+1) <= Vmin2
                        Vol_OP_Real_2(t+1) = Vmin2;
                        IP2 = 1;
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
                    %%%%% This fixes the volume of the dams between what we
                    % define as maximum and minimum during the simulation
                    % of the HJB equation (UP).
                    
                    % Here we redefine the costs to no considere the
                    % penalization when we compute the accumulated one.
                    dd_OP = dd;
                    d_OP = @(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t);
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t);
                    
                    % Accumulated cost over the optimal path.
                    instantTotalCost(t) = dt*TMax*(Controls_OP(t,:)*...
                        [d_OP(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),uV1,uV2,uV3,uV4,uA,t,Prec_Cont(t,:));...
                        olCV2(uV2,t);olCV3(uV3,t)] + Controls_OP(t,1:12)*dd_OP*Controls_OP(t,1:12)' +...
                        H0_OP(uV1,uV2,uV3,uV4,t,Prec_Cont(t,:)));
                    Optimal_Cost(t+1) = Optimal_Cost(t) + instantTotalCost(t);

                    Battery_Cost = Battery_Cost + dt*TMax*Controls_OP(t,7)*olCA(uA);

                else
                    
                    % Here we redefine the costs to no considere the
                    % penalization when we compute the accumulated one.
                    dd_OP = dd;
                    d_OP = @(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t,Conts) d(V1,V2,V3,uV1,uV2,uV3,uV4,uA,t);
                    H0_OP = @(uV1,uV2,uV3,uV4,t,Conts) H0(uV1,uV2,uV3,uV4,t);
                    instantTotalCost(t) = dt*TMax*(Controls_OP(t,:)*...
                    [d_OP(Vol_OP_Real_1(t),Vol_OP_Real_2(t),Vol_OP_Real_3(t),uV1,uV2,uV3,uV4,uA,t,Prec_Cont(t,:));...
                    olCV2(uV2,t);olCV3(uV3,t)] + Controls_OP(t,1:12)*dd_OP*Controls_OP(t,1:12)' +...
                    H0_OP(uV1,uV2,uV3,uV4,t,Prec_Cont(t,:)));
                    
                end
                
                if WhatToDo ~= 7
                    Power_H1_OP(t) = Power_H1(Controls_OP(t,1),Controls_OP(t,2),Vol_OP_Real_1(t),Vol_OP_Real_2(t));
                    Power_H2_OP(t) = Power_H2(Controls_OP(t,3),Controls_OP(t,4),Vol_OP_Real_2(t),Vol_OP_Real_3(t));
                    Power_H3_OP(t) = Power_H3(Controls_OP(t,5),Vol_OP_Real_3(t),t);
                    Power_H4_OP(t) = Power_H4(Controls_OP(t,6),Vol_OP_Real_4(t),t);
                    Power_Battery_OP(t) = Power_Battery(Controls_OP(t,7));
                    Power_F1_OP(t) = MaxTermicas(1)*(Controls_OP(t,8));
                    Power_F2_OP(t) = MaxTermicas(2)*(Controls_OP(t,9));
                    Power_F3_OP(t) = MaxTermicas(3)*(Controls_OP(t,10));
                    Power_F4_OP(t) = MaxTermicas(4)*(Controls_OP(t,11));
                    Power_F5_OP(t) = MaxTermicas(5)*(Controls_OP(t,12));
                elseif WhatToDo == 7
                    Power_H1_OP = histPowerBon;
                    Power_H2_OP = histPowerBay;
                    Power_H3_OP = histPowerPal;
                    Power_H4_OP = histPowerSG;
                    Power_Battery_OP = zeros(1,length(time));
                    Power_F1_OP = histPowerMot;
                    Power_F2_OP = histPowerPTA;
                    Power_F3_OP = histPowerPTB;
                    Power_F4_OP = histPowerCTR;
                    Power_F5_OP = zeros(1,length(time));
                end
                
                waterBoneteOverTime(t)      = FlowMax1(Vol_OP_Real_1(t),Vol_OP_Real_2(t))*Controls_OP(t,1) + ...
                    MaxSpill1(Vol_OP_Real_1(t),Vol_OP_Real_2(t))*Controls_OP(t,2);
                waterBaygorriaOverTime(t)   = FlowMax2(Vol_OP_Real_2(t),Vol_OP_Real_3(t))*Controls_OP(t,3) + ...
                    MaxSpill2(Vol_OP_Real_2(t),Vol_OP_Real_3(t))*Controls_OP(t,4);
                waterPalmarOverTime(t)      = FlowMax3(Vol_OP_Real_3(t),t)*Controls_OP(t,5);
                waterSGOverTime(t)          = FlowMax4*Controls_OP(t,6);
                virtualBaygorriaOverTime(t) = MaxV2*Controls_OP(t,13);
                virtualPalmarOverTime(t)    = MaxV3*Controls_OP(t,14);
                
                % We multiply by dt*TMax which is the number of seconds per discretization.
                % Then WaterBonete and WaterVirtualBay are in m^3.
                WaterBonete(t)           = dt*TMax * waterBoneteOverTime(t); 
                WaterBaygorria(t)        = dt*TMax * waterBaygorriaOverTime(t); 
                WaterPalmar(t)           = dt*TMax * waterPalmarOverTime(t); 
                WaterSG(t)               = dt*TMax * waterSGOverTime(t); 
                WaterVirtualBaygorria(t) = dt*TMax * virtualBaygorriaOverTime(t);
                WaterVirtualPalmar(t)    = dt*TMax * virtualPalmarOverTime(t);

            end
            
            for t = 1:length(time)-1
                futureCostOP(t) = sum(instantTotalCost(t:end-1));
                futureRealCostOP(t) = sum(realInstantCost(t:end-1));
                instantCostHJB(t) = uPlot(t) - uPlot(t+1);
            end
            
            if WhatToDo ~= 4 && WhatToDo ~= 6 && WhatToDo ~= 7
                for i = 1:length(time)
                    Spot(i) = lambda{i}.ineqnonlin*MWhTokW; % In USD/kW.
                    spotTimesPH1(i) = Spot(i)*dPH1overFlowPlot(i); % In USD/(m^3/s).
                    spotTimesPH2(i) = Spot(i)*dPH2overFlowPlot(i); % In USD/(m^3/s).
                    spotTimesPH1(i) = spotTimesPH1(i) * dt*TMax;
                    spotTimesPH2(i) = spotTimesPH2(i) * dt*TMax;
                end
            end    

            % ==================== 'Others' Section ====================>>

            if ComputePlots == 1

                Format_Long = 1;
                
                normBalance = 1000*3000;
                if Format_Long == 0
                    set(0,'CurrentFigure',101); clf(101); set(gcf,'Position',[10 10 900 420],'name','Power Balance');
                else
                	set(0,'CurrentFigure',101); clf(101); set(gcf,'Position',[10 10 750*2 300*2],'name','Power Balance');
                end
                grid minor; hold on;
                if UseBattery == 1
                    hh = area(time,[BioMasPower;SolarPower;WindPower;Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                    Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP;Power_Battery_OP.*(Power_Battery_OP > 0)]'/normBalance);
                    hh(13).FaceColor = [1 1 0]; % Battery.
                else
                    hh = area(time,[BioMasPower;SolarPower;WindPower;Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                    Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP]'/normBalance);
                end
                hh(1).FaceColor = [1 .5 0]; % Biomasa.
                hh(2).FaceColor = [1 1 .5]; % Solar.
                hh(3).FaceColor = [0 1 0]; % Wind.
                hh(4).FaceColor = [.5 .5 1]; % SG.
                hh(5).FaceColor = [.5 .8 1]; % Bonete.
                hh(6).FaceColor = [.2 .2 1]; % Baygorria.
                hh(7).FaceColor = [.8 .8 1]; % Palmar.
                hh(8).FaceColor = [1 .6 .6]; % T1.
                hh(9).FaceColor = [1 .3 .5]; % T2.
                hh(10).FaceColor = [145/255 0 211/255;]; % T3.
                hh(11).FaceColor = [1 0 0]; % T4.
                hh(12).FaceColor = [0 0 0]; % T5.
                P = plot(time,D/normBalance,'k'); P.LineWidth = 3;
                if UseBattery == 1
                    P = plot(time,(D-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0))/normBalance,'c');
                    P.LineWidth = 3;
                end
                P = plot(time,(D+Export-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0))/normBalance,'m');
                P.LineWidth = 3;
                if UseBattery == 1
                    if Format_Long == 0
                      %  legend('Biomass Power','Solar Power','Wind Power','Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)',...
                      %      'Motores Batlle (Fossil)','PTA (Fossil)','PTB (Fossil)','CTR (Fossil)','Failure (Fossil)','Battery','Demand','Demand + Battery','Dem. + Battery + Exp.','location','eastoutside');
                        legend('Biomass Power','Solar Power','Wind Power','Salto Grande','Bonete','Baygorria','Palmar',...
                            'Motores Batlle','PTA','PTB','CTR','Failure','Battery','Demand','Demand + Battery','Dem. + Battery + Exp.','location','eastoutside');
                    else
                       % legend('Biomass Power','Solar Power','Wind Power','Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)',...
                        %    'Motores Batlle (Fossil)','PTA (Fossil)','PTB (Fossil)','CTR (Fossil)','Failure (Fossil)','Battery','Demand','Demand + Battery','Dem. + Battery + Exp.','location','eastoutside');
                    	legend('Biomass Power','Solar Power','Wind Power','Salto Grande','Bonete','Baygorria','Palmar',...
                            'Motores Batlle','PTA','PTB','CTR','Failure','Battery','Demand','Demand + Battery','Dem. + Battery + Exp.','location','eastoutside');
                    end
                else
   %                 legend('Biomass Power','Solar Power','Wind Power','Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)',...
    %                    'Motores Batlle (Fossil)','PTA (Fossil)','PTB (Fossil)','CTR (Fossil)','Failure (Fossil)','Demand','Demand + Exp.','location','eastoutside');
                    legend('Biomass Power','Solar Power','Wind Power','Salto Grande','Bonete','Baygorria','Palmar',...
                        'Motores Batlle','PTA','PTB','CTR','Failure','Demand','Demand + Exp.','location','eastoutside');
                end
                xlabel('Time');
                title('Simulated Daily Power Balance')
                %ylabel('Power (MW)');
                ylabel('Power');
                xlim([0 max(time)]);
                set(findall(gcf,'-property','FontSize'),'FontSize',24);
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
                P = plot(time,Optimal_Cost/NormCost); P.LineWidth = 1;
                P = plot(time,futureCostOP/NormCost); P.LineWidth = 1;
                P = plot(time,futureRealCostOP/NormCost); P.LineWidth = 1;
                P = plot(time,uPlot/NormCost); P.LineWidth = 1;
                P = plot(time,ones(1,length(time))*Final_Cost/NormCost); P.LineWidth = 1;
                P = plot(time,accumRealInstantCost/NormCost); P.LineWidth = 1;
                line([0,1],[0,0],'Color','red','LineStyle','--')
                grid on; grid minor;
                legend('Optimal Path Cost','Optimal Path Future Cost','Optimal Path Real Future Cost','HJB Cost Optimal Path',...
                    'HJB Final Cost','Real cost Optial Path',...
                    'Zero reference','location','southeast');
                xlabel('Time'); ylabel('Hundred Thousand Dollars');
                title('Accumulated Simulated Cost');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',103); clf(103); set(gcf,'name','Historic Dams Controls');
                hold on; grid on; grid minor;
                title('Historical Dams Controls');
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
                title('Historical Fossil Fuel Stations Controls');
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
                plot(time,Vol_OP_Real_1);
                ylim([Vmin1 Vmax1]);
                grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Bonete over Historical Path');
                hold on;
                xlim([0 max(time)]);
                for i = 2:length(V1)-1
                    line([0,1],[V1(i),V1(i)],'Color','green','LineStyle','--')
                end
                legend('Volume over Historical Path');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',107); clf(107); set(gcf,'name','Volume Palmar');
                plot(time,Vol_OP_Real_3);
                ylim([Vmin3 Vmax3]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Palmar over Historical Path');
                xlim([0 max(time)]);
                for i = 2:length(V3)-1
                    line([0,1],[V3(i),V3(i)],'Color','green','LineStyle','--')
                end
                legend('Volume over Historical Path');
                set(gca,'Box','on');

                set(0,'CurrentFigure',108); clf(108); set(gcf,'name','Volume Salto Grande');
                plot(time,Vol_OP_Real_4);
                ylim([Vmin4 Vmax4]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Salto Grande over Historical Path');
                xlim([0 max(time)]);
                for i = 2:length(V4)-1
                    line([0,1],[V4(i),V4(i)],'Color','green','LineStyle','--')
                end
                legend('Volume over Historical Path');
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

                set(0,'CurrentFigure',110); clf(110); set(gcf,'name','Costs in the Optimal Path');
                hold on;
                P = plot(time,CostControls(1,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(3,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(5,:)*MWhTokW); P.LineWidth = 1;
                P = plot(time,CostControls(6,:)*MWhTokW); P.LineWidth = 1;
                if UseBattery == 1
                    P = plot(time,CostControls(7,:)*MWhTokW); P.LineWidth = 1;
                end
                grid on; grid minor;
                title('Costs over Historical Path');
                if UseBattery == 1
                legend('Turbine flow Bonete','Turbine flow Baygorria','Turbine flow Palmar','Turbine flow Salto Grande','Battery','location','southeast');
                else
                legend('Turbine flow Bonete','Turbine flow Baygorria','Turbine flow Palmar','Turbine flow Salto Grande','location','southeast');
                end
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
                if UseBattery == 1
                P = plot(time,CostControls(7,:)*MWhTokW);
                P.LineWidth = 1;
                end
                grid on; grid minor;
                title('Costs over Historical Path (Fuel)');
                if UseBattery == 1
                legend(Termicas{1:4},'Battery'); xlabel('Time'); ylabel('USD/MWh');
                else
                legend(Termicas{1:4}); xlabel('Time'); ylabel('USD/MWh');
                end
                ylim([0.9*min([min(CostControls(8,:)),min(CostControls(9,:)),min(CostControls(10,:)),min(CostControls(11,:))])...
                    1.1*max([max(CostControls(8,:)),max(CostControls(9,:)),max(CostControls(10,:)),max(CostControls(11,:))])]*MWhTokW);
                set(gca,'Box','on');

                set(0,'CurrentFigure',113); clf(113); set(gcf,'name','Cost of using Spillage');
                hold on;
                P = plot(time,CostControls(2,:)); P.LineWidth = 1;
                P = plot(time,CostControls(4,:)); P.LineWidth = 1;
                grid on; grid minor;
                title('Cost of using Spillage');
                legend('Bonete Spillage','Baygorria Spillage');
                xlabel('Time'); ylabel('USD/hm^3');
                set(gca,'Box','on');

                set(0,'CurrentFigure',114); clf(114); set(gcf,'name','Instant Derivatives Values');
                hold on;
                plot(time,uV1plot/NormCost); plot(time,uV2plot/NormCost); plot(time,uV3plot/NormCost);
                plot(time,uV4plot/NormCost); plot(time,uAplot/NormCost);
                grid on; grid minor;
                Leg2 = legend('$\frac{\partial u}{\partial \hat{V}^{(1)}}$','$\frac{\partial u}{\partial \hat{V}^{(2)}}$','$\frac{\partial u}{\partial \hat{V}^{(3)}}$',...
                    '$\frac{\partial u}{\partial \hat{V}^{(4)}}$','$\frac{\partial u}{\partial A}$','location','southeast');
                set(Leg2,'interpreter', 'latex','FontSize',14);
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
                if UseBattery == 1
                X = [sum(BioMasPower),sum(SolarPower),sum(WindPower),sum(Power_H1_OP),sum(Power_H2_OP),...
                    sum(Power_H3_OP),sum(Power_H4_OP),sum(Power_F1_OP),sum(Power_F2_OP),sum(Power_F3_OP),...
                    sum(Power_F4_OP),sum(Power_F5_OP),sum(Power_Battery_OP.*(Power_Battery_OP > 0))];
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
                else
                X = [sum(BioMasPower),sum(SolarPower),sum(WindPower),sum(Power_H1_OP),sum(Power_H2_OP),...
                    sum(Power_H3_OP),sum(Power_H4_OP),sum(Power_F1_OP),sum(Power_F2_OP),sum(Power_F3_OP),...
                    sum(Power_F4_OP),sum(Power_F5_OP)];
                txt = {'Biomass: ','Solar Power: ','Wind Power: ',[Dams{1},': '],[Dams{2},': '],...
                    [Dams{3},': '],[Dams{4},': '],[Termicas{1},': '],[Termicas{2},': '],...
                    [Termicas{3},': '],[Termicas{4},': '],[Termicas{5},': ']};
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
                0 0 0];
                end
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
                delete(ax.Children(1:2:end));
                combinedtxt = strcat(txt,percentValues');
                legend(combinedtxt,'Location','eastoutside');
                colormap(colors);
                set(findall(gcf,'-property','FontSize'),'FontSize',16);
                set(gca,'Box','on');

                set(0,'CurrentFigure',117); clf(117); set(gcf,'name','Costs Distribution');
                X = [sum(costControllable(1,:)),sum(costControllable(2,:)),sum(costControllable(3,:)),sum(costControllable(4,:)),...
                    sum(costControllable(5,:)),sum(costControllable(6,:)),sum(costControllable(7,:)),sum(costControllable(8,:)),sum(costControllable(10,:))];
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
                title('Controllable Historical Costs Distribution');
                pText = findobj(p2,'Type','text');
                percentValues = get(pText,'String');
                ax = gca;
                delete(ax.Children(1:2:end));
                combinedtxt = strcat(txt,percentValues');
                legend(combinedtxt,'Location','eastoutside');
                colormap(colors);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',123); clf(123); set(gcf,'name','Costs Histogram');
                X = [sum(costControllable(1,:)),sum(costControllable(2,:)),sum(costControllable(3,:)),sum(costControllable(4,:)),...
                    sum(costControllable(5,:)),sum(costControllable(6,:)),sum(costControllable(7,:)),sum(costControllable(8,:)),...
                    sum(costControllable(10,:))]/NormCost;
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
                title('Historical Controllable Costs Distribution');
                legend(txt,'Location','eastoutside');
                ylabel('Hundred Thousand Dollars');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',124); clf(124); set(gcf,'Position',[10 10 900 420],'name','Cost Optimal Path');
                grid minor; hold on;
                hh = area(time,[accumContCost(1,:);accumContCost(2,:);accumContCost(3,:);accumContCost(4,:);accumContCost(5,:);...
                    accumContCost(6,:);accumContCost(7,:);accumContCost(8,:);accumContCost(9,:);accumContCost(10,:);accumContCost(11,:)]'/NormCost);
                hh(1).FaceColor = [.5 .5 1]; % Bonete.
                hh(2).FaceColor = [.5 .8 1]; % Baygorria.
                hh(3).FaceColor = [.2 .2 1]; % Palmar.
                hh(4).FaceColor = [.8 .8 1]; % SG.
                hh(5).FaceColor = [1 .6 .6]; % T1.
                hh(6).FaceColor = [1 .3 .5]; % T2.
                hh(7).FaceColor = [145/255 0 211/255;]; % T3.
                hh(8).FaceColor = [1 0 0]; % T4.
                hh(9).FaceColor = [0 0 0]; % T5.
                hh(10).FaceColor = [0 1 0.5]; % V2.
                hh(11).FaceColor = [0.5 1 0]; % V3.
                txt = [Dams(:)',Termicas(:)','Virtual 2','Virtual 3'];
                legend(txt,'location','eastoutside');
                xlabel('Time');
                title('Historical Controllable Accumulated Cost')
                ylabel('Hundred Thousand Dollars');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',125); clf(125); set(gcf,'Position',[10 10 900 420],'name','Cost Optimal Path');
                grid minor; hold on;
                hh = area(time,[accumCost(1,:);accumCost(2,:);accumCost(3,:);accumCost(4,:);accumCost(5,:);...
                    accumCost(6,:);accumCost(7,:);accumCost(8,:);accumCost(9,:);accumCost(10,:);accumCost(11,:)]'/NormCost,'FaceAlpha',0.4);
                hh(1).FaceColor = [.5 .5 1]; % Bonete.
                hh(2).FaceColor = [.5 .8 1]; % Baygorria.
                hh(3).FaceColor = [.2 .2 1]; % Palmar.
                hh(4).FaceColor = [.8 .8 1]; % SG.
                hh(5).FaceColor = [1 .6 .6]; % T1.
                hh(6).FaceColor = [1 .3 .5]; % T2.
                hh(7).FaceColor = [145/255 0 211/255;]; % T3.
                hh(8).FaceColor = [1 0 0]; % T4.
                hh(9).FaceColor = [0 0 0]; % T5.
                hh(10).FaceColor = [0 1 0.5]; % V2.
                hh(11).FaceColor = [0.5 1 0]; % V3.
                P = plot(time,Optimal_Cost/NormCost); P.LineWidth = 2;
                if WhatToDo == 7 || WhatToDo == 4 || WhatToDo == 3 || WhatToDo == 2 || WhatToDo == 6
                    P = plot(time,ones(1,length(time))*Final_Cost); P.LineWidth = 2;
                else
                	P = plot(time,ones(1,length(time))*Final_Cost/NormCost); P.LineWidth = 2;
                end
                legend(txt{:},'Optimal Path Cost','HJB Final Cost','location','eastoutside');
                xlabel('Time');
                title('Accumulated Cost')
                ylabel('Hundred Thousand Dollars');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',126); clf(126); set(gcf,'name','Costs Histogram');
                X = [sum(totalInstantCost(1,:)),sum(totalInstantCost(2,:)),sum(totalInstantCost(3,:)),sum(totalInstantCost(4,:)),...
                    sum(costControllable(5,:)),sum(costControllable(6,:)),sum(costControllable(7,:)),sum(costControllable(8,:)),sum(costControllable(10,:))]/NormCost;
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
                title('Historical Costs Distribution');
                legend(txt,'Location','eastoutside');
                ylabel('Hundred Thousand Dollars');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',127); clf(127); set(gcf,'name','Volume Baygorria');
                plot(time,Vol_OP_Real_2);
                ylim([Vmin2 Vmax2]);
                grid on; grid minor;
                xlabel('Time');
                ylabel('Volume');
                title('Volume Baygorria over Historical Path');
                xlim([0 max(time)]);
                for i = 2:length(V2)-1
                    line([0,1],[V2(i),V2(i)],'Color','green','LineStyle','--')
                end
                legend('Volume over Historical Path');
                set(gca,'Box','on');
                
                if WhatToDo ~= 4 && WhatToDo ~= 6 && WhatToDo ~= 7
                    set(0,'CurrentFigure',118); clf(118); set(gcf,'name','Spot Price');
                    plot(time,Spot);grid minor;
                    legend('Spot Price');
                    title('Spot Price');
                    xlabel('Time'); ylabel('USD/MWh');
                    set(gca,'Box','on');
                end

                if Format_Long == 0
                    set(0,'CurrentFigure',120); clf(120); set(gcf,'Position',[10 10 900 420],'name','Power Balance Controllable');
                else
                	set(0,'CurrentFigure',120); clf(120);% set(gcf,'Position',[10 10 750 750],'name','Power Balance Controllable');
                end
                grid minor; hold on;
                if UseBattery == 1
                    hh = area(time,[Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                        Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP;Power_Battery_OP.*(Power_Battery_OP > 0)]'/normBalance);
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
                else
                    hh = area(time,[Power_H4_OP;Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_F1_OP;...
                    Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_F5_OP]'/normBalance);
                    hh(1).FaceColor = [.5 .5 1]; % Bonete.
                    hh(2).FaceColor = [.5 .8 1]; % Baygorria.
                    hh(3).FaceColor = [.2 .2 1]; % Palmar.
                    hh(4).FaceColor = [.8 .8 1]; % SG.
                    hh(5).FaceColor = [1 .6 .6]; % T1.
                    hh(6).FaceColor = [1 .3 .5]; % T2.
                    hh(7).FaceColor = [145/255 0 211/255;]; % T3.
                    hh(8).FaceColor = [1 0 0]; % T4.
                    hh(9).FaceColor = [0 0 0]; % T5.
                end
                if Format_Long == 0
                    if UseBattery == 1
                        legend(Dams{4},Dams{1:3},Termicas{1:5},'Battery','location','eastoutside');
                    else
                    	legend(Dams{4},Dams{1:3},Termicas{1:5},'location','eastoutside');
                    end
                else
                    if UseBattery == 1
                        P = plot(time,(D+Export-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0)-BioMasPower-SolarPower-WindPower)/normBalance,'m');
                        P.LineWidth = 2;
                        legend(Dams{4},Dams{1:3},Termicas{1:5},'Battery','Dem. + Battery + Exp.','location','SouthEast');
                    else
                        P = plot(time,(D+Export-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0)-BioMasPower-SolarPower-WindPower)/normBalance,'m');
                        P.LineWidth = 2;
                    	legend(Dams{4},Dams{1:3},Termicas{1:5},'Demand + Exp.','location','SouthEast');
                    end
                end
                xlabel('Time');
                title('Simulated Daily Power Balance (Controllable)')
                ylabel('Power');
                xlim([0 max(time)]);
                set(findall(gcf,'-property','FontSize'),'FontSize',14);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',121); clf(121); set(gcf,'name','Costs Distribution');
                X = [sum(costControllable(1,:))+sum(costControllable(2,:))+sum(costControllable(3,:))+sum(costControllable(4,:)),...
                    sum(costControllable(5,:))+sum(costControllable(6,:))+sum(costControllable(7,:))+sum(costControllable(8,:))+sum(costControllable(10,:))];
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
                title('Hydro-Fuel Cost Ratio');
                pText = findobj(p2,'Type','text');
                percentValues = get(pText,'String');
                ax = gca;
                delete(ax.Children(1:2:end));
                combinedtxt = strcat(txt,percentValues');
                legend(combinedtxt,'Location','southoutside');
                colormap(colors);
                set(findall(gcf,'-property','FontSize'),'FontSize',20);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',122); clf(122); set(gcf,'name','Virtual controls over the Optimal Path');
                hold on; grid minor;
                title('Controls over the Optimal Path');
                xlabel('Time');
                ylabel('Control');
                P = plot(time,Controls_OP(:,13)); P.LineWidth = 1;
                P = plot(time,Controls_OP(:,14)); P.LineWidth = 1;
                legend('Virtual Control 1','Virtual Control 2');
                xlim([0 max(time)]); ylim([0 1]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',128); clf(128); set(gcf,'name','Battery Capacity');
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
                saveas(gcf,[pwd '/',['Simulations/Simulation_',name],'/Extra_',num2str(i)],'epsc');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',129); clf(129); set(gcf,'name','Instant Cost');
                hold on;
                plot(time,instantTotalCost/NormCost,'LineWidth',2);
                plot(time,instantCostHJB/NormCost);
                plot(time,realInstantCost/NormCost);
                grid minor;
                legend('Instant Cost OP','Instant Cost HJB','Instant Real Cost OP','location','southeast');
                xlabel('Time'); ylabel('Hundred Thousand Dollars');
                title('Instant Cost');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',130); clf(130); set(gcf,'name','Marginal Costs');
                hold on;
                plot(time,ones(1,length(time))*CH1*1e6); % The 1e6 passes from USD/m^3 to USD/hm^3.
                plot(time,spotTimesPH1*1e6); % The 1e6 passes from USD/m^3 to USD/hm^3.
                grid minor;
                legend('Instant Cost OP','Instant Cost HJB','Instant Real Cost OP','location','southeast');
                xlabel('Time'); ylabel('USD/hm$^3$','Interpreter','latex');
                title('Instant Cost');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',131); clf(131); set(gcf,'name','Costs in the Optimal Path');
                hold on;
                P = plot(time,CV2InstantCost); P.LineWidth = 1;
                P = plot(time,CV3InstantCost); P.LineWidth = 1;
                grid on; grid minor;
                title('Costs over Historical Path (Virtual Controls)');
                legend('Instant Cost Virtual Baygorria','Instant Cost Virtual Palmar')
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',132); clf(132); set(gcf,'name','Water Exchange');
                hold on;
                plot(time,-WaterBonete,'LineWidth',1);
                plot(time,-WaterBaygorria,'LineWidth',1);
                plot(time,-WaterPalmar,'LineWidth',1);
                plot(time,-WaterSG,'LineWidth',1);
                plot(time,WaterVirtualBaygorria,'LineWidth',1);
                plot(time,WaterVirtualPalmar,'LineWidth',1);
                plot(time,ones(1,length(time))*IT1(1)*dt*TMax,'LineWidth',1);
                plot(time,ones(1,length(time))*IT2(1)*dt*TMax,'LineWidth',1);
                plot(time,ones(1,length(time))*IT3(1)*dt*TMax,'LineWidth',1);
                plot(time,ones(1,length(time))*IT4(1)*dt*TMax,'LineWidth',1);
                grid minor;
                title('Water Exchange');
                ylabel('m^3');
                legend('Used Bonete','Used Baygorria','Used Palmar','Used SG','Virtual Baygorria',...
                    'Virtual Palmar', 'Inflow Bonete', 'Inflow Baygorria', 'Inflow Palmar', 'Inflow SG');
                set(gca,'Box','on');
                
                set(0,'CurrentFigure',133); clf(133); set(gcf,'name','Water Balance');
                hold on;
                plot(time,ones(1,length(time))*IT1(1)*dt*TMax - WaterBonete,'LineWidth',1);
                plot(time,ones(1,length(time))*IT2(1)*dt*TMax - WaterBaygorria + WaterVirtualBaygorria,'LineWidth',1);
                plot(time,ones(1,length(time))*IT3(1)*dt*TMax - WaterPalmar + WaterVirtualPalmar,'LineWidth',1);
                plot(time,ones(1,length(time))*IT4(1)*dt*TMax - WaterSG,'LineWidth',1);
                grid minor;
                title('Water Balance');
                ylabel('m^3');
                legend('Balance Bonete','Balance Baygorria','Balance Palmar','Balance SG');
                set(gca,'Box','on');

                if SavePlots == 1
                    for i = Plot_12:Plot_22
                        set(0,'CurrentFigure',i);
                        if WhatToDo == 7
                            saveas(gcf,[pwd '/',['Simulations/Real_',name],'/',num2str(i)],'epsc');
                        else
                            saveas(gcf,[pwd '/',['Simulations/Simulation_',name],'/',num2str(i)],'epsc');
                        end
                    end
                    i = 11;
                    set(0,'CurrentFigure',i);
                    if WhatToDo == 7
                        saveas(gcf,[pwd '/',['Simulations/Real_',name],'/',num2str(i)],'epsc');
                    else
                        saveas(gcf,[pwd '/',['Simulations/Simulation_',name],'/',num2str(i)],'epsc');
                    end
                end

            end

            Final_Cost = accumRealInstantCost(end-1);
            if GradNum == 1000 
                save([pwd '/','Simulations/Simulation_Controls_',name,'_Day_',num2str(Day),'.mat'],'Controls_OP');
            end
            toc

            if WhatToDo == 6

                for i = 0:lenLambda21-1
                    if i == lenLambda21-1
                        gradientBaygorria(i+1) = sum(WaterVirtualBaygorria(Z21+i*D21:T-1) - ...
                            WaterBonete(Z21+i*D21-Z21+1:T-1-Z21+1));
                        D21Last = (T-1-Z21+1) - (Z21+i*D21-Z21+1) + 1;
                    else
                        gradientBaygorria(i+1) = sum(WaterVirtualBaygorria(Z21+i*D21:Z21+(i+1)*D21-1) - ...
                            WaterBonete(1+i*D21:(i+1)*D21));
                    end
                end
                for i = 0:lenLambda32-1
                    if i == lenLambda32-1
                        gradientPalmar(i+1) = sum(WaterVirtualPalmar(Z32+i*D32:T-1) - ...
                            WaterBaygorria(Z32+i*D32-Z32+1:T-1-Z32+1));
                        D32Last = (T-1-Z32+1) - (Z32+i*D32-Z32+1) + 1;
                    else
                        gradientPalmar(i+1) = sum(WaterVirtualPalmar(Z32+i*D32:Z32+(i+1)*D32-1) - ...
                            WaterBaygorria(1+i*D32:(i+1)*D32));
                    end
                end
               
                figure(1000);
                set(gcf,'name','Lagrangian Multipliers');
                lamTime = linspace(0,1+t32,length([Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]));
                plot(lamTime,[Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]*1e6,'b--o');
                hold on;
                plot(lamTime,Lambda32*1e6,'r--o');
                grid minor;
                title('Approximated Lagrangian Multipliers');
                ylabel('USD/hm$^3$','Interpreter','latex');
                xlabel('Time');
                legend('$\hat{\lambda}_{2,1}(t)$','$\hat{\lambda}_{3,2}(t)$','interpreter', 'latex');
                xlim([0 max(lamTime)]);
                                
                figure(1001); 
                set(gcf,'name','Subgradient');
                plot(linspace(0.25,1,length(gradientBaygorria)),gradientBaygorria/1e6,'b--o');
                hold on; 
                title('Subgradient Bonete-Baygorria');
                xlabel('Time');
                ylabel('hm^3');
                xlim([0.25 max(time)]);
                legend('Subgradient $\lambda_{2,1}(t)$','interpreter', 'latex');
                set(gca,'Box','on');
                grid minor;
                
                figure(1002); 
                set(gcf,'name','Subgradient');
                plot(linspace(10/24,1,length(gradientPalmar)),gradientPalmar/1e6,'b--o');
                hold on; 
                title('Subgradient Baygorria-Palmar');
                xlabel('Time');
                ylabel('hm^3');
                xlim([10/24 max(time)]);
                legend('Subgradient $\lambda_{3,2}(t)$','interpreter', 'latex');
                set(gca,'Box','on');
                grid minor;
                
                figure(1003); 
                set(gcf,'name','Real Water');
                plot(time,waterBoneteOverTime,'b--o');
                hold on;
                plot(time,virtualBaygorriaOverTime,'r--o');
                title('Real Bonete - Virtual Baygorria');
                legend('Total Bonete Flow','Virtual Baygorria Flow');
                xlabel('Time');
                ylabel('m^3/s');
                xlim([0 max(time)]);
                
                set(gca,'Box','on');
                grid minor;
                
                figure(1004); 
                set(gcf,'name','Real Water');
                plot(time,waterBaygorriaOverTime,'b--o');
                hold on;
                plot(time,virtualPalmarOverTime,'r--o');
                title('Real Baygorria - Virtual Palmar');
                legend('Total Baygorria Flow','Virtual Palmar Flow');
                xlabel('Time');
                ylabel('m^3/s');
                xlim([0 max(time)]);
                set(gca,'Box','on');
                grid minor;
                
                Final_Cost = [accumRealInstantCost(end-1),gradientBaygorria*1e-6,gradientPalmar*1e-6];

%                 for i = 1000:1004
%                     set(0,'CurrentFigure',i);
%                     saveas(gcf,[pwd '/',['Simulations/Simulation_',name],'/',num2str(i)],'epsc');
%                 end
                
                disp(['D21 = ',num2str(D21),' and TMax*dt*D21 = ',num2str(D21*dt*TMax)]);
                disp(['D32 = ',num2str(D32),' and TMax*dt*D32 = ',num2str(D32*dt*TMax)]);
                disp(['D21Last = ',num2str(D21Last),' and TMax*dt*D21Last = ',num2str(D21Last*dt*TMax)]);
                disp(['D32Last = ',num2str(D32Last),' and TMax*dt*D32Last = ',num2str(D32Last*dt*TMax)]);
                
            end
        
        end
        
    end
    
    %% End.

end