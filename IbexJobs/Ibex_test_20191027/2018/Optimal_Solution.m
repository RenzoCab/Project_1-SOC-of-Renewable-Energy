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
    
    % NEW!!!:
    if 0 == exist([pwd '/Simulations/Simulation_',name],'dir')
        mkdir([pwd '/Simulations/Simulation_',name]);
    end
    if 0 == exist([pwd '/Simulations/Simulation_OPT_',name],'dir')
        mkdir([pwd '/Simulations/Simulation_OPT_',name]);
    end

    load('MaxFlowCoeff.mat');
    load('Day_20180302.mat');
    
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

    NT = 2^(Expan_Time); % Discretizations in time.
    NV1 = 2^(Expan-2); % Discretizations of Bonete.
    NV2 = 2^(Expan); % Discretizations of Baygorria.
    NV3 = 2^(Expan-2); % Discretizations of Palmar.
    NV4 = 2^(Expan+1); % Discretizations of Salto Grande.
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

    % NEW!!!: I added still the extra +0.1 in Bay. And the else includes
    % WhatToDo == 6.
    if (WhatToDo == 1) % Solve HJB.
        Khe1 = Matrix{3}(3); % In USD/MWh. Cost per energy of Bonete.
        Khe2 = Matrix{3}(4)+0.1; % In USD/MWh. Cost per energy of Baygorria.
        Khe3 = Matrix{3}(2); % In USD/MWh. Cost per energy of Palmar.
        Khe4 = Matrix{3}(1); % In USD/MWh. Cost per energy of SG.
        water1zero = 1; water2zero = 1; water3zero = 1; water4zero = 1;
    else % Else is only the last post processing. No WhatToDo == 6.
        minWater1 = 2e-2; minWater2 = 1e-2; minWater3 = 1.5e-2; minWater4 = 0.5e-2;
        Khe1 = max(Matrix{3}(3),minWater1); % In USD/MWh. Cost per energy of Bonete.
        Khe2 = max(Matrix{3}(4)+0.1,minWater2); % In USD/MWh. Cost per energy of Baygorria.
        Khe3 = max(Matrix{3}(2),minWater3); % In USD/MWh. Cost per energy of Palmar.
        Khe4 = max(Matrix{3}(1),minWater4); % In USD/MWh. Cost per energy of SG.
        water1zero = (Matrix{3}(3) >= minWater1);
        water2zero = (Matrix{3}(4) >= minWater2);
        water3zero = (Matrix{3}(2) >= minWater3);
        water4zero = (Matrix{3}(1) >= minWater4);
    end
    
    if artificialWaterValue
%         Khe1 = 0; Khe2 = 0; Khe3 = 0; Khe4 = 0; 
%         Khe1 = 1;
%         Khe2 = 0;
%         Khe3 = 0;
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
    
	tic

        for t = length(time):-1:1

            % ==================== Controls Section ====================>>
            disp(num2str(100*t/length(time)));
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
                
                % NEW!!!:
                Controls_Aux(index,:) = [Quad_FMC_Lambda_8(Q(V1(ind_V1),V2(ind_V2),V3(ind_V3),t),b(V1(ind_V1),V2(ind_V2),V3(ind_V3),...
                        V4(ind_V4),t),c(t),d(V1(ind_V1),V2(ind_V2),V3(ind_V3),uV1_PF(index),v2_aux,v3_aux,uV4_PF(index),...
                        Aux_uA,t),dd,options,WS,Aux_Lims,InfTermicas,name);Aux_V2;Aux_V3];

            end

            Controls_T = reshape(Controls_Aux,LV1,LV2,LV3,LV4,LA,14);
        
        end
        
        toc
        
end