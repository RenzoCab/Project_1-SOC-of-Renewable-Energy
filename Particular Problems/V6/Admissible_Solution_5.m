function [HJB_Initial_Cost] = Admissible_Solution_5(name,WhatToDo,WarmStart,UseBattery,...
    ComputePlots,NewPlots,SavePlots,Expan,Expan_Time,UseLamb21,UseLamb32,...
    Lambda21,Lambda32)

    Phi_V2 = 0.2; % Value of the virtual control for t <= Z21.
    Phi_V3 = 0.2; % Value of the virtual control for t <= Z32.

    Derivatives = [1,1,1,1]; % For debugging.

    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',1e4,...
        'FunctionTolerance',1e-10,'MaxIterations',1e3,'MaxFunEvals',1e7,'Algorithm','sqp',...
        'ConstraintTolerance',1e-5,'StepTolerance',1e-15);

    if NewPlots == 1
        close all;
        NeededPlots(1,14,0);
        if WhatToDo == 2 || WhatToDo == 3
            NeededPlots(101,118,0);
        end
    elseif NewPlots ~= 0
        disp('Choose NewPlots 0 or 1.');
        return;
    end
    
    MWhTokW = 3.6e6; % From MWh to kW (energy in 1 h to power per second).
    m3Tohm3 = 1e6;  % From m^3 to hm^3.
    Termicas = {'Central Batlle','Punta del Tigre: 1 to 6','Punta del Tigre: 7 to 8',...
        'Central Térmica Respaldo'};
    Dams = {'Bonete','Baygorria','Palmar','Salto Grande'};
    MaxTermicas = [60,288,48,208]*1000; % kW.
    MinTermicas = [6,90,0.6,40]*1000*0; % kW. We multiply by 0 to make all them 0.
    CostTermicas = [82,86,88,103]; % USD/MWh.
    CostTermicasNorm = CostTermicas/MWhTokW; % USD/kW.
    InfTermicas = MinTermicas./MaxTermicas;

    NT = 2^Expan_Time; % Discretizations in time.
    NV1 = 2^Expan; % Discretizations of Bonete.
    NV2 = 2^Expan; % Discretizations of Baygorria.
    NV3 = 2^Expan; % Discretizations of Palmar.
    NV4 = 2^Expan; % Discretizations of Salto Grande.
    NA = 2^(Expan+3); % Discretizations of Salto Grande.

    % ==================== Time ====================>>>

    absTime = 0:1/24:1; % Is needed for interpolation.
    Tmax = 1; % Simulation final time.
    Tmin = 0; % Simulation initial time.
    dt = (Tmax-Tmin)/NT;
    time = Tmin:dt:Tmax; % Simulation time.

    % ==================== Parameters (Demand, Water and Fuel) ====================>>

    H1 = @(V1) (-3.74)*(V1)^2 + (16.7)*(V1) + 67.7; % Height of the water as a function of the normal volume (Bonete).
    % H2 = @(V2) (-3.25)*(V2).^2 + (12.1)*(V2) + 45.6; % Height of the water as a function of the normal volume (Baygorria).
    H2 = 54; % We fix the level of Baygorria to its nominal level.
    H3 = @(V3) (-3.76)*(V3)^2 + (16.9)*(V3) + 26.8; % Height of the water as a function of the normal volume (Palmar).
    H4 = @(V4) (-19.8)*(V4)^2 + (51.5)*(V4) + 3.79; % Height of the water (m) as a function of the normal volume (Salto Grande).
    h04 = 5.1; % Level after Salto Grande. Now is a constant, but may depend of the time.
    h03 = 15; % Level after Palmar. Now is a constant, but may depend of the time.

    RealMaxFlow1 = 640; % In m^3/s.
    RealMaxFlow2 = 828; % In m^3/s.
    RealMaxFlow3 = 1373; % In m^3/s.
    RealMaxFlow4 = 4200; % In m^3/s.

%     FlowMax1 = @(V1) (H1(V1)-H2 <= 21)*(12.4*(H1(V1)-H2) + 447) + (H1(V1)-H2 > 21)*(-11.2*(H1(V1)-H2) + 941); % In m^3/s.
%     FlowMax2 = @(V3) (H2-H3(V3) <= 14)*(37.5*(H2-H3(V3)) + 376) + (H2-H3(V3) > 14)*(-81.3*(H2-H3(V3)) + 2039); % In m^3/s.
%     FlowMax3 = @(V3) (H3(V3)-h03 <= 21)*(46.5*(H3(V3)-h03) + 638) + (H3(V3)-h03 > 21)*(-68.1*(H3(V3)-h03) + 3272); % In m^3/s.
    FlowMax1 = @(V1) 0.8*(H1(V1)-H2 <= 21)*(7.6739*(H1(V1)-H2) + 518) + (H1(V1)-H2 > 21)*(-7.8004*(H1(V1)-H2) + 843)*0.8; % In m^3/s.
    FlowMax2 = @(V3) (H2-H3(V3) <= 14)*(37.4862*(H2-H3(V3)) + 303) + (H2-H3(V3) > 14)*(-81.3*(H2-H3(V3)) + 1966); % In m^3/s.
    FlowMax3 = @(V3) (H3(V3)-h03 <= 21)*(46.4714*(H3(V3)-h03) + 304) + (H3(V3)-h03 > 21)*(-68.1*(H3(V3)-h03) + 2938); % In m^3/s.
    FlowMax4 = 4200; % In m^3/s.

    MaxSpill1 = @(V1) 3060 - FlowMax1(V1); % Otherwise is 2420 m^3/s.
    MaxSpill2 = @(V3) 3474 - FlowMax2(V3); % Otherwise is 2646 m^3/s.
    MaxSpill3 = @(V3) 4160 - FlowMax3(V3); %Otherwise is 2787 m^3/s.
    MaxSpill4 = 12600; % In m^3/s.

    MaxRel1 = 3060; % In m^3/s (640+2420).
    MaxRel2 = 3474; % In m^3/s (828+2646).
    MaxRel3 = 4160; % In m^3/s (1373+2787).
    MaxRel4 = FlowMax4 + MaxSpill4; % In m^3/s.

    MaxV2 = MaxRel1;
    MaxV3 = MaxRel2;

    % Here we use an unreal natural inflow equal to half of the maximum turbine
    % flow, for all time. Also, we compute the maximum and minimum posible
    % volumes of the dams given that natural inputs.
    IT1 = @(t) RealMaxFlow1/5; % In m^3/s.
    IT2 = @(t) RealMaxFlow2/20; % In m^3/s.
    IT3 = @(t) RealMaxFlow3/20; % In m^3/s.
    IT4 = @(t) RealMaxFlow4/5; % In m^3/s.

    VolMax1 = 9.5e9; % In m^3.
    VolMax2 = 617e6; % In m^3.
    VolMax3 = 2823e6; % In m^3.
    VolMax4 = 5e9; % In m^3.

    VolMin1 = 1925e6; % In m^3.
    VolMin2 = 424e6; % In m^3.
    VolMin3 = 1767e6; % In m^3.
    VolMin4 = 3470e6; % In m^3.

    TMax = 24*3600; % In s. Number of seconds in a day.

    f1 = @(t,V1,T,S) (TMax/VolMax1)*(IT1(t)-FlowMax1(V1)*T-MaxSpill1(V1)*S);
    f3 = @(t,V3,T,S,V) (TMax/VolMax3)*(IT3(t)-FlowMax3(V3)*T-MaxSpill3(V3)*S+MaxV3*V);
    f4 = @(t,V4,T,S) (TMax/VolMax4)*(IT4(t)-FlowMax4*T-MaxSpill4*S);

    FuelMax = 2e6; % In kW.
    DMax = 2e6; % In kW. The maximum possible value of the demand.

    Vini1 = 0.8; % Initial condition between 0.2 and 1.
    Vini2 = 0.8; % Initial condition between 0.7 and 1.
    Vini3 = 0.8; % Initial condition between 0.5 and 1.
    Vini4 = 0.9; % Initial condition between 0.6 and 1.

    t21 = 1/4; % Delay between Bonete and Baygorria.
    Z21 = floor(length(time)*t21); % Discrete delay between Bonete and Baygorria.
    t32 = 1/6; % Delay between Baygorria and Palmar.
    Z32 = floor(length(time)*t32); % Discrete delay between Baygorria and Palmar.
    t31 = 10/24; % Delay between Bonete and Palmar.
    Z31 = floor(length(time)*t31); % Discrete delay between Bonete and Palmar.

    if UseLamb21 == 0 || (UseLamb21 == 1 && length(Lambda21) == 1)
        Lambda21 = zeros(1,length(time)+Z21);
    end
    if UseLamb32 == 0 || (UseLamb32 == 1 && length(Lambda32) == 1)
        Lambda32 = zeros(1,length(time)+Z32);
    end
    Lambda21 = interp1(0:1/(length(Lambda21)-1):1,Lambda21,0:1/(length(time)+Z21-1):1);
    Lambda32 = interp1(0:1/(length(Lambda32)-1):1,Lambda32,0:1/(length(time)+Z32-1):1);
    if ComputePlots == 1
        set(0,'CurrentFigure',11);
        clf(11);
        if t21 > t32
            lamTime = 0:dt:1+t21;
        else
            lamTime = 0:dt:1+t31;
        end
        plot(lamTime,[Lambda21,zeros(1,length(Lambda32)-length(Lambda21))]*1e6);
        hold on; grid on;
        plot(lamTime,[Lambda32,zeros(1,length(Lambda21)-length(Lambda32))]*1e6);
        title('Approximated Lagrangian multipliers');
        ylabel('USD/hm$^3$','Interpreter','latex');
        xlabel('Normalized Time');
        Leg = legend('$\hat{\lambda}_{2,1}(t)$','$\hat{\lambda}_{3,2}(t)$');
        set(Leg, 'interpreter', 'latex');
        xlim([0 max(lamTime)]);
    elseif ComputePlots ~= 0
        disp('Choose ComputePlots between 0 or 1.');
        return;
    end

    sigH1 = 0; % Sigma from the SDE of Bonete.
    sigH2 = 0; % Sigma from the SDE of Baygorria.
    sigH3 = 0; % Sigma from the SDE of Palmar.
    sigH4 = 0; % Sigma from the SDE of Salto Grande.
    NormCost = 1e5; % The normalization is over 100.000 dollars.

    Vmax1 = min((Vini1*VolMax1+(RealMaxFlow1)*TMax),VolMax1)/VolMax1;
    Vmin1 = max((Vini1*VolMax1-(RealMaxFlow1)*TMax),VolMin1)/VolMax1;
    Vmax2 = min((Vini2*VolMax2+(MaxRel2)*TMax),VolMax2)/VolMax2;
    Vmin2 = max((Vini2*VolMax2-(MaxRel2)*TMax),VolMin2)/VolMax2;
    Vmax3 = min((Vini3*VolMax3+(RealMaxFlow3)*TMax),VolMax3)/VolMax3;
    Vmin3 = max((Vini3*VolMax3-(RealMaxFlow3)*TMax),VolMin3)/VolMax3;
    Vmax4 = min((Vini4*VolMax4+(FlowMax4)*TMax),VolMax4)/VolMax4;
    Vmin4 = max((Vini4*VolMax4-(FlowMax4)*TMax),VolMin4)/VolMax4;

    load('D01012017.mat'); % Demand from 01/01/2017. It loads a 25x1 vector named D01012017.

    dV1 = (Vmax1-Vmin1)/NV1;
    dV2 = (Vmax2-Vmin2)/NV2;
    dV3 = (Vmax3-Vmin3)/NV3;
    dV4 = (Vmax4-Vmin4)/NV4;

    V1 = Vmin1:dV1:Vmax1; % Discretized volume of Bonete.
    V2 = Vmin2:dV2:Vmax2; % Discretized volume of Baygorria.
    V3 = Vmin3:dV3:Vmax3; % Discretized volume of Palmar.
    V4 = Vmin4:dV4:Vmax4; % Discretized volume of Salto Grande.

    [vplot1_4,vplot3_4] = meshgrid(V1,V3); % For plotting surfaces.
    [vplot1_3,vplot4_3] = meshgrid(V1,V4); % For plotting surfaces.
    [vplot3_1,vplot4_1] = meshgrid(V3,V4); % For plotting surfaces.

    D = D01012017'*1000*1.2;
    D = interp1(absTime,D,time)+DMax*0.15*sin(4*pi*time);

    % ==================== Plot Demand (plot 6) ====================>>

    if ComputePlots == 1
        clf(1);
        set(0,'CurrentFigure',1);
        hold on;
        xlabel('Normalized Time');
        title('Normalized Demand');
        grid on;
        plot(time,D/DMax);
        ylim([0 1]);
        pause(0.001);
    end

    % ==================== Battery ====================>>

    AMax = 1.4e5; % In kWh.
    PAMax = 1e5; % In kW.
    K_A = (PAMax/AMax)/3600; % Normalizing constant.
    sigA = 0; % Sigma from the SDE of the battery.
    A0 = 0.5;
    NAMax = 1;
    NAMin = 0;
    olPn = @(A) min(1,(A-NAMin)/(K_A*dt*TMax));
    ulPn = @(A) -min(0.35,(NAMax-A)/(K_A*dt*TMax));

%     Par = [0,1,0.3,20];
%     olPn = @(A) (Par(1)*exp(Par(3)*Par(4))+Par(2).*exp(Par(4).*A)) ./ (exp(Par(3)*Par(4))+exp(Par(4).*A));
%     Par2 = [0,1,0.8,20];
%     ulPn = @(A) ((Par2(1)*exp(Par2(3)*Par2(4))+Par2(2).*exp(Par2(4).*A)) ./ (exp(Par2(3)*Par2(4))+exp(Par2(4).*A)) - 1)*0.35;

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
        clf(13);
        set(0,'CurrentFigure',13);
        hold on; grid on; grid minor;
        xlabel('A');
        ylabel('Control');
        title('Battery Control Vs Charge');
        grid on;
        plot(A,olPn(A));
        plot(A,ulPn(A));
        legend('Maximum Control','Minimum Control');
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
    ControlsFMC = zeros(length(V1),length(V3),length(V4),length(A),13); % T1,S1,T2,S2,T3,T4,A,F1,F2,F3,F4,V2,V3.
    LV1 = length(V1); % Length of V1.
    LV3 = length(V3); % Length of V3.
    LV4 = length(V4); % Length of V4.
    LA = length(A); % Length of A.
    TL = LV1*LV3*LV4*LA; % Total number of elements in the matrix.
    Controls_T = ControlsFMC;
    UV = {};

    % ==================== More parameters (Water and Fuel) ====================>>

    cT1 = 1.02;
    cS1 = 3.9;
    d1 = @(tur,spill) cT1*tur+cS1*spill; % Variation in the high due of the released flow (Bonete).

    cT2 = 0.57;
    cS2 = 1.8;
    d2 = @(tur,spill) cT2*tur+cS2*spill; % Variation in the high due of the released flow (Baygorria).

    cT3 = 3.77;
    cS3 = 7.65;
    d3 = @(tur,spill) cT3*tur+cS3*spill; % Variation in the high due of the released flow (Palmar).

    cT4 = 5.34;
    cS4 = cT4*MaxSpill4/FlowMax4;
    d4 = @(tur,spill) cT4*tur+cS4*spill; % Variation in the high due of the released flow (Salto Grande).

    eta = 8.72; % Efficiency.

    Khe1 = 2.0; % In USD/MWh. Cost per energy of Bonete.
    Khe2 = 0; % In USD/MWh.
    Khe3 = 2.3; % In USD/MWh.
    Khe4 = 4.3; % In USD/MWh.
    Kf = 13; % In USD/MWh. Cost per energy of the fuel.

    CH1 = Khe1*eta*(H1(Vini1)-d1(1,0)-H2)/(MWhTokW); % Water's value of Bonete.
    CH2 = Khe2*eta*(H2-d2(1,0)-H3(Vini3))/(MWhTokW); % Water's value of Baygorria.
    CH3 = Khe3*eta*(H3(Vini3)-d3(1,0)-h03)/(MWhTokW); % Water's value of Palmar.
    CH4 = Khe4*eta*(H4(Vini4)-d4(1,0)-h04)/(MWhTokW); % Water's value of Salto Grande.
    CF = Kf/(MWhTokW); % USD/kW.

    % We divide over 3600 to pass from hours to seconds and over 1000 to move
    % from Mega to Kilo.

%     olCF = CF*FuelMax;
    olCT1 = @(uV1,t,V1) FlowMax1(V1)*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCS1 = @(uV1,t,V1) MaxSpill1(V1)*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCT2 = @(uV2,t,V3) FlowMax2(V3)*(CH2 - Lambda32(t+Z32) - uV2/VolMax2 + Lambda21(t));
    olCS2 = @(uV2,t,V3) MaxSpill2(V3)*(CH2 - Lambda32(t+Z32) - uV2/VolMax2 + Lambda21(t));
    olCT3 = @(uV3,V3) FlowMax3(V3)*(CH3 - uV3/VolMax3);
    olCT4 = @(uV4) FlowMax4*(CH4 - uV4/VolMax4);
    olCA = @(uA) -K_A*uA -KA*1000000*0;
    H0 = @(uV1,uV2,uV3,uV4,t) IT1(t)*uV1/VolMax1 + IT3(t)*uV3/VolMax3 + IT4(t)*uV4/VolMax4 - Lambda21(t)*IT2(t);
    olCV3 = @(uV3,t) MaxRel3*(Lambda32(t) + uV3/VolMax3);

    k1 = @(t) IT2(t)/MaxV2;
    k2 = @(V3) FlowMax2(V3)/MaxV2;
    k3 = @(V3) MaxSpill2(V3)/MaxV2;
    k = @(t,V3) [k1(t),k2(V3),k3(V3)];

    % Hydraulic power: PH = T^2*K2 + T*K1 + T*S*K3.

    K11 = @(V1) eta*FlowMax1(V1)*(H1(V1)-H2);
    K12 = @(V3) eta*FlowMax2(V3)*(H2-H3(V3));
    K13 = @(V3) eta*FlowMax3(V3)*(H3(V3)-h03);
    K14 = @(V4) eta*FlowMax4*(H4(V4)-h04);

    K21 = @(V1) -eta*FlowMax1(V1)*cT1;
    K22 = @(V3) -eta*FlowMax2(V3)*cT2;
    K23 = @(V3) -eta*FlowMax3(V3)*cT3;
    K24 = -eta*FlowMax4*cT4;

    K31 = @(V1) -eta*FlowMax1(V1)*cS1;
    K32 = @(V3) -eta*FlowMax2(V3)*cS2;
    K33 = @(V3) -eta*FlowMax3(V3)*cS3;
    K34 = -eta*FlowMax4*cS4;

    d = @(V1,V3,uV1,uV2,uV3,uV4,uA,t) [olCT1(uV1,t,V1),olCS1(uV1,t,V1),...
        olCT2(uV2,t,V3),olCS2(uV2,t,V3),olCT3(uV3,V3),olCT4(uV4),olCA(uA),...
        CostTermicasNorm.*MaxTermicas]';
    dd = diag([0,0,0,0,0,0,KA*1000000*0.000,0,0,0,0]);
    rho = KA*1000000*0.07;
    b = @(V1,V3,V4) [K11(V1),0,K12(V3),0,K13(V3),K14(V4),PAMax,MaxTermicas]';
    Q = @(V1,V3) [K21(V1) K31(V1)/2 0 0 0 0 0 0 0 0 0;
        K31(V1)/2 0 0 0 0 0 0 0 0 0 0;
        0 0 K22(V3) K32(V3)/2 0 0 0 0 0 0 0;
        0 0 K32(V3)/2 0 0 0 0 0 0 0 0;
        0 0 0 0 K23(V3) 0 0 0 0 0 0;
        0 0 0 0 0 K24 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0]; % Is fixed.

    c = @(t) -D(t);

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

                                if (IT3(t+1)-FlowMax3(V3(v3))*ControlsFMC_Aux(v1,v3,v4,a,5)+MaxV3*ControlsFMC_Aux(v1,v3,v4,a,13) >= 0) && Derivatives(2)
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
                                    if IT3(t+1)-FlowMax3(V3(v3))*ControlsFMC_Aux(v1,v3,v4,a,5)+MaxV3*ControlsFMC_Aux(v1,v3,v4,a,13) >= 0
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

%                                 f1 = @(t,V1,T,S) (TMax/VolMax1)*(IT1(t)-FlowMax1(V1)*T-MaxSpill1(V1)*S);
%                                 f3 = @(t,V3,T,S,V) (TMax/VolMax3)*(IT3(t)-FlowMax3(V3)*T-MaxSpill3(V3)*S+MaxV3*V);
%                                 f4 = @(t,V4,T,S) (TMax/VolMax4)*(IT4(t)-FlowMax4*T-MaxSpill4*S);
%                                 fA = @(C) TMax*K_A*C;

                            CFL(v1,v3,v4,a) = 1 / (abs(f1(t,V1(v1),ControlsFMC_Aux(v1,v3,v4,a,1),ControlsFMC_Aux(v1,v3,v4,a,2)))/dV1 + ...
                                abs(f3(t,V3(v3),ControlsFMC_Aux(v1,v3,v4,a,5),0,ControlsFMC_Aux(v1,v3,v4,a,13)))/dV3 + ...
                                abs(f4(t,V4(v4),ControlsFMC_Aux(v1,v3,v4,a,6),0))/dV4 + ...
                                abs(fA(ControlsFMC_Aux(v1,v3,v4,a,7)))/dA);
                            CFL_v1(v1,v3,v4,a) = abs(f1(t,V1(v1),ControlsFMC_Aux(v1,v3,v4,a,1),ControlsFMC_Aux(v1,v3,v4,a,2)))/dV1;
                            CFL_v3(v1,v3,v4,a) = abs(f3(t,V3(v3),ControlsFMC_Aux(v1,v3,v4,a,5),0,ControlsFMC_Aux(v1,v3,v4,a,13)))/dV3;
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
            Controls_Aux = zeros(TL,13);
            uA_PF = uA(:);
            Warm_Start = reshape(Controls_T,LV1,LV3,LV4,LA,13);
            Warm_Start = Warm_Start(:,:,:,:,1:11);

            if WarmStart == 0
                Warm_Start(Warm_Start ~= 0) = 0;
            elseif WarmStart ~= 1
                disp('Choose WarmStart 0 or 1.');
                return;
            end

            parfor index = 1:TL

                if (UseLamb32 == 0) || (olCV3(uV3(index),t) >= 0 && t > Z32)
                    Aux_V3 = 0;
                elseif t > Z32
                    Aux_V3 = 1;
                else
                    Aux_V3 = Phi_V3;
                end

                if UseLamb21 == 1
                    if t > Z21
                        Aux_V2 = 0;
                        PhiV2Fixed = 0;
                    else
                        Aux_V2 = Phi_V2;
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

                Controls_Aux(index,:) = [Quad_FMC_Lambda_5(Q(V1(ind_V1),V3(ind_V3)),b(V1(ind_V1),V3(ind_V3),...
                        V4(ind_V4)),c(t),d(V1(ind_V1),V3(ind_V3),uV1_PF(index),0,uV3_PF(index),uV4_PF(index),...
                        Aux_uA,t),dd,k(t,V3(ind_V3)),options,PhiV2Fixed,Aux_V2,WS,Aux_Lims,InfTermicas);Aux_V3];

            end

            Controls_T = reshape(Controls_Aux,LV1,LV3,LV4,LA,13);

            % ==================== Time Step ====================>>

            if t ~= length(time)
                for v1 = 1:length(V1)
                    for v3 = 1:length(V3)
                        for v4 = 1:length(V4)
                            for a = 1:length(A)

                                Controls(1:13) = ControlsFMC_Aux(v1,v3,v4,a,:);

                                HAM(v1,v3,v4,a) = TMax*(Controls*[d(V1(v1),V3(v3),uV1(v1,v3,v4,a),...
                                    0,uV3(v1,v3,v4,a),uV4(v1,v3,v4,a),uA(v1,v3,v4,a),t);0;...
                                    olCV3(uV3(v1,v3,v4,a),t)]+H0(uV1(v1,v3,v4,a),0,uV3(v1,v3,v4,a),uV4(v1,v3,v4,a),t));

                                u(v1,v3,v4,a) = u(v1,v3,v4,a) + dt*HAM(v1,v3,v4,a);

                            end
                        end
                    end
                end
            end

            % ==================== Plotting over time ====================>>

            if ComputePlots == 1

                set(0,'CurrentFigure',1);
                vline(time(t),'r',' ');
                title(['Normalized Demand (time=',num2str(t),')']);
                if SavePlots == 1 && t == 1
                    saveas(gcf,['Demand_',num2str(Expan)],'epsc');
                end

                set(0,'CurrentFigure',2);
                clf(2); hold on;
                for index = 1:floor(length(V4)/2):length(V4)
                    surf(vplot1_4,vplot3_4,squeeze(u(:,:,index,1)/NormCost)');
                    title(['Optimal Cost (empty Battery) at time = ',num2str(t)]);
                end
                xlabel('Vol. Bonete'); ylabel('Vol. Palmar');
                xlim([Vmin1 Vmax1]); ylim([Vmin3 Vmax3]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',3);
                clf(3); hold on;
                for index = 1:floor(length(V3)/2):length(V3)
                    surf(vplot1_3,vplot4_3,squeeze(u(:,index,:,1)/NormCost)');
                    title(['Optimal Cost (empty Battery) at time = ',num2str(t)]);
                end
                xlabel('Vol. Bonete'); ylabel('Vol. Salto Grande');
                xlim([Vmin1 Vmax1]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',4);
                clf(4);
                hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                    surf(vplot3_1,vplot4_1,squeeze(u(index,:,:,1)/NormCost)');
                    title(['Optimal Cost (empty Battery) at time = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',5);
                clf(5); hold on;
                for index = 1:floor(length(V4)/2):length(V4)
                        surf(vplot1_4,vplot3_4,squeeze(HAM(:,:,index,1)/NormCost)');
                        title(['Hamiltonian (empty Battery) at time = ',num2str(t)]);
                end
                xlabel('Vol. Bonete'); ylabel('Vol. Palmar');
                xlim([Vmin1 Vmax1]); ylim([Vmin3 Vmax3]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',6);
                clf(6); hold on;
                for index = 1:floor(length(V3)/2):length(V3)
                        surf(vplot1_3,vplot4_3,squeeze(HAM(:,index,:,1)/NormCost)');
                        title(['Hamiltonian (empty Battery) at time = ',num2str(t)]);
                end
                xlabel('Vol. Bonete'); ylabel('Vol. Salto Grande');
                xlim([Vmin1 Vmax1]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',7);
                clf(7); hold on;
                for index = 1:floor(length(V1)/2):length(V1)
                        surf(vplot3_1,vplot4_1,squeeze(HAM(index,:,:,1)/NormCost)');
                        title(['Hamiltonian (empty Battery) at time = ',num2str(t)]);
                end
                xlabel('Vol. Palmar'); ylabel('Vol. Salto Grande');
                xlim([Vmin3 Vmax3]); ylim([Vmin4 Vmax4]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',8);
                clf(8); hold on;
                for index = 1:floor(length(V4)/2):length(V4)
                    surf(vplot3_A,vplotA_3,squeeze(uV4(index,:,index,:)/NormCost));
                end
                xlabel('Vol. Palmar'); ylabel('Charge Battery');
                title(['u_{V4} at time = ',num2str(t)]);
                xlim([Vmin3 Vmax3]); ylim([NAMin NAMax]);
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

                set(0,'CurrentFigure',12);
                clf(12); hold on;
                for i = 1:floor(length(V4)/2):length(V4)
                    surf(vplot4_A,vplotA_4,squeeze(u(i,i,:,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['Optimal Cost at time = ',num2str(t)]);
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

                set(0,'CurrentFigure',14);
                clf(14); hold on;
                for index = 1:length(V4)
                    surf(vplot4_A,vplotA_4,squeeze(uA(index,:,index,:)/NormCost));
                end
                xlabel('Vol. Salto Grande'); ylabel('Charge Battery');
                title(['u_{A} at time = ',num2str(t)]);
                xlim([Vmin4 Vmax4]); ylim([NAMin NAMax]);
                zlim auto; box; view(140,20);

            end

            fprintf('Discretizations = %d. Completed: %.2f.\n',NV4,100*(1-t/length(time)));
            pause(0.1);

        end

        IP1 = min(abs(V1-Vini1)) == abs(V1-Vini1);
        IP3 = min(abs(V3-Vini3)) == abs(V3-Vini3);
        IP4 = min(abs(V4-Vini4)) == abs(V4-Vini4);
        IPA = min(abs(A-A0)) == abs(A-A0);
        HJB_Initial_Cost = u(IP1,IP3,IP4,IPA);

        ToSave = {HJB_Initial_Cost,UV,CFL_Plot};
        save(['Simulation_',name,'.mat'],'ToSave');

        toc
        
    end
    
    %% Optimal Path Section:
    
    if WhatToDo == 2 || WhatToDo == 3

        tic

        load(['Simulation_',name,'.mat']);
        HJB_Initial_Cost = ToSave{1};
        UV = ToSave{2};
        CFL_Plot = ToSave{3};
        
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
        Cost = zeros(9,length(time));
        uV1plot = zeros(1,length(time));
        uV3plot = zeros(1,length(time));
        uV4plot = zeros(1,length(time));
        CFL_plot = zeros(1,length(time));
        Controls_OP = zeros(length(time),10);
        CostControls = zeros(11,length(time));
        Battery_OP = zeros(1,length(time));
        Battery_Real = zeros(1,length(time));
        Power_Battery_OP = zeros(1,length(time));
        Prev_SG = zeros(1,length(time)+1);
        Prev_Batt = zeros(1,length(time)+1);
        Prev_F4 = zeros(1,length(time)+1);
        Power_F1_OP = zeros(1,length(time));
        Power_F2_OP = zeros(1,length(time));
        Power_F3_OP = zeros(1,length(time));
        Power_F4_OP = zeros(1,length(time));
        Extra_Cost = 0;
        
        uAplot = zeros(1,length(time));
        Battery_OP = A(IPA);
        Battery_Real(1) = A(IPA);
        x0 = zeros(11,1);
        
        Power_H1 = @(x1,x2,V1) eta*x1*FlowMax1(V1)*(H1(V1)-d1(x1,x2)-H2); % In kW.
        Power_H2 = @(x3,x4,V3) eta*x3*FlowMax2(V3)*(H2-d2(x3,x4)-H3(V3)); % In kW.
        Power_H3 = @(x5,V3) eta*x5*FlowMax3(V3)*(H3(V3)-d3(x5,0)-h03); % In kW.
        Power_H4 = @(x6,V4) eta*x6*FlowMax4*(H4(V4)-d4(x6,0)-h04); % In kW.
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
                    Real2 = Phi_V2;
                else
                    Real2 = 0;
                end
            end
            
            if t == 1
                dd_OP = dd;
                d_OP = @(V1,V3,uV1,uV2,uV3,uV4,uA,t,ContSG,ContA,ContF4) d(V1,V3,uV1,uV2,uV3,uV4,uA,t);
                H0_OP = @(uV1,uV2,uV3,uV4,t,ContSG,ContA,ContF4) H0(uV1,uV2,uV3,uV4,t);
            else
                dd_OP = dd + diag([0,0,0,0,0,rho/dt,rho/dt,0,0,0,rho/dt]);
                d_OP = @(V1,V3,uV1,uV2,uV3,uV4,uA,t,ContSG,ContA,ContF4) d(V1,V3,uV1,uV2,uV3,uV4,uA,t) + ...
                [0,0,0,0,0,-2*rho*ContSG/dt,-2*rho*ContA/dt,0,0,0,-2*rho*ContF4/dt]';
                H0_OP = @(uV1,uV2,uV3,uV4,t,ContSG,ContA,ContF4) H0(uV1,uV2,uV3,uV4,t) + (rho*ContA^2)/dt + (rho*ContF4^2)/dt + (rho*ContSG^2)/dt;
                Added_Cost = @(New_ContSG,New_ContA,New_ContF4,ContSG,ContA,ContF4) (New_ContSG^2+New_ContA^2+New_ContF4^2)*rho/dt - ...
                    2*rho/dt*(ContSG*New_ContSG+ContA*New_ContA+ContF4*New_ContF4) + (ContSG^2+ContA^2+ContF4^2)*rho/dt;
            end
            
            Controls_OP(t,1:12) = Quad_FMC_Lambda_5(Q(Vol_OP_Real_1(t),Vol_OP_Real_3(t)),...
            	b(Vol_OP_Real_1(t),Vol_OP_Real_3(t),Vol_OP_Real_4(t)),c(t),...
            	d_OP(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t,Prev_SG(t),Prev_Batt(t),Prev_F4(t)),dd_OP,k(t,Vol_OP_Real_3(t)),...
                options,1,Real2,x0,ContLims,InfTermicas);

            Phi_V2_Aux(t) = (Controls_OP(t,1)*FlowMax1(Vol_OP_Real_1(t)) + Controls_OP(t,2)*MaxSpill1(Vol_OP_Real_1(t))) / MaxV2;
            Phi_V3_Aux(t) = (Controls_OP(t,3)*FlowMax2(Vol_OP_Real_3(t)) + Controls_OP(t,4)*MaxSpill2(Vol_OP_Real_3(t))) / MaxV3;
            x0(1:11) = Controls_OP(t,1:11);
            Prev_Batt(t+1) = Controls_OP(t,7);
            Prev_F4(t+1) = Controls_OP(t,11);
            Prev_SG(t+1) = Controls_OP(t,6);
            
            if UseLamb32 == 1
                if t > Z32
                    Controls_OP(t,13) = Phi_V3_Aux(t-Z32);
                else
                    Controls_OP(t,13) = Phi_V3;
                end
            else
                Controls_OP(t,13) = 0;
            end
            
            Aux_Cost = Controls_OP(t,1:11).*d(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t)';
            
            Cost(1,t) = Aux_Cost(1)+Aux_Cost(2);
            Cost(2,t) = Aux_Cost(3)+Aux_Cost(4);
            Cost(3,t) = Aux_Cost(5);
            Cost(4,t) = Aux_Cost(6);
            Cost(5,t) = Aux_Cost(8);
            Cost(6,t) = Aux_Cost(9);
            Cost(7,t) = Aux_Cost(10);
            Cost(8,t) = Aux_Cost(11);
            Cost(9,t) = Aux_Cost(7);
            
            CostControls(:,t) = d(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t)./...
               [Power_H1(1,Controls_OP(t,2),Vol_OP_Real_1(t)),MaxSpill1(Vol_OP_Real_1(t))*dt*TMax/m3Tohm3,...
               Power_H2(1,Controls_OP(t,4),Vol_OP_Real_3(t)),MaxSpill2(Vol_OP_Real_3(t))*dt*TMax/m3Tohm3,...
               Power_H3(1,Vol_OP_Real_3(t)),Power_H4(1,Vol_OP_Real_4(t)),PAMax,MaxTermicas]';
           % This cost is in USD/kW for all but the spillages. For the
           % spillages is USD/Control.
           
            % ==================== Dynamics and Powers Section ====================>>

            if t ~= length(time)
                                    
                Vol_OP_Real_1(t+1) = Vol_OP_Real_1(t) + (dt*TMax/VolMax1)*(IT1(t)-...
                    Controls_OP(t,1)*FlowMax1(Vol_OP_1(t))-Controls_OP(t,2)*MaxSpill1(Vol_OP_1(t)));
                Vol_OP_Real_3(t+1) = Vol_OP_Real_3(t) + (dt*TMax/VolMax3)*(IT3(t)-...
                    Controls_OP(t,5)*FlowMax3(Vol_OP_3(t))+Controls_OP(t,13)*MaxV3);
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
                if Battery_Real(t+1) >= 1
                    Battery_Real(t+1) = 1;
                    IPA = length(A)-1;
                elseif Battery_Real(t+1) < 0
                    Battery_Real(t+1) = 0;
                    IPA = 1;
                end
                
                if t == 1
                    Optimal_Cost(t+1) = Optimal_Cost(t) + dt*TMax*(Controls_OP(t,:)*...
                        [d(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t);0;olCV3(uV3,t)] +...
                        Controls_OP(t,1:11)*dd*Controls_OP(t,1:11)' + H0_OP(uV1,uV2,uV3,uV4,t,Prev_SG(t),Prev_Batt(t),Prev_F4(t)));
                else
                    Optimal_Cost(t+1) = Optimal_Cost(t) + dt*TMax*(Controls_OP(t,:)*...
                    [d_OP(Vol_OP_Real_1(t),Vol_OP_Real_3(t),uV1,0,uV3,uV4,uA,t,Prev_SG(t),Prev_Batt(t),Prev_F4(t));0;olCV3(uV3,t)] +...
                    Controls_OP(t,1:11)*dd_OP*Controls_OP(t,1:11)' + H0_OP(uV1,uV2,uV3,uV4,t,Prev_SG(t),Prev_Batt(t),Prev_F4(t)) -...
                    Added_Cost(Controls_OP(t,6),Controls_OP(t,7),Controls_OP(t,11),Prev_SG(t),Prev_Batt(t),Prev_F4(t)));
                    Extra_Cost = Extra_Cost + dt*TMax*Added_Cost(Controls_OP(t,6),Controls_OP(t,7),Controls_OP(t,11),Prev_SG(t),Prev_Batt(t),Prev_F4(t))/NormCost;
                end
                
            end

            % ==================== 'Others' Section ====================>>

            
            Power_H1_OP(t) = Power_H1(Controls_OP(t,1),Controls_OP(t,2),Vol_OP_Real_1(t));
            Power_H2_OP(t) = Power_H2(Controls_OP(t,3),Controls_OP(t,4),Vol_OP_Real_3(t));
            Power_H3_OP(t) = Power_H3(Controls_OP(t,5),Vol_OP_Real_3(t));
            Power_H4_OP(t) = Power_H4(Controls_OP(t,6),Vol_OP_Real_4(t));
            Power_Battery_OP(t) = Power_Battery(Controls_OP(t,7));
            Power_F1_OP(t) = MaxTermicas(1)*(Controls_OP(t,8));
            Power_F2_OP(t) = MaxTermicas(2)*(Controls_OP(t,9));
            Power_F3_OP(t) = MaxTermicas(3)*(Controls_OP(t,10));
            Power_F4_OP(t) = MaxTermicas(4)*(Controls_OP(t,11));
                       
        end

        if ComputePlots == 1

            set(0,'CurrentFigure',101); clf(101);
            grid on; grid minor; hold on;
            hh = area(time,[Power_H1_OP;Power_H2_OP;Power_H3_OP;Power_H4_OP;Power_F1_OP;...
                Power_F2_OP;Power_F3_OP;Power_F4_OP;Power_Battery_OP.*(Power_Battery_OP > 0)]');
            P = plot(time,D,'k'); P.LineWidth = 1;
            P = plot(time,D-squeeze(Controls_OP(:,7))'*PAMax.*(squeeze(Controls_OP(:,7))' < 0),'g');
            P.LineWidth = 1;
            legend(Dams{1:4},Termicas{1:4},'Battery','Demand','Demand + Battery','location','southeast');
            hh(1).FaceColor = [.5 .5 1];
            hh(2).FaceColor = [.5 .8 1];
            hh(3).FaceColor = [.2 .2 1];
            hh(4).FaceColor = [.8 .8 1];
            hh(5).FaceColor = [1 .6 .6];
            hh(6).FaceColor = [1 .3 .5];
            hh(7).FaceColor = [1 .5 .3];
            hh(8).FaceColor = [1 0 0];
            hh(9).FaceColor = [1 1 .5];
            xlabel('Normalized Time');
            title('Power Balance')
            ylabel('Power (kW)');
            xlim([0 max(time)]);

            set(0,'CurrentFigure',102); clf(102);
            hold on;
            plot(time,Optimal_Cost/NormCost);
            plot(time,ones(1,length(time))*HJB_Initial_Cost/NormCost);
            grid on; grid minor;
            legend('Accumulated Cost','HJB Final Cost','location','southeast');
            xlabel('Time'); ylabel('Cost');
            title('Costs Comparison');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',103); clf(103);
            hold on; grid on; grid minor;
            title('Controls over the Time');
            xlabel('Time');
            ylabel('Control');
            for q = 1:6
                P = plot(time,Controls_OP(:,q));
                P.LineWidth = 1;
            end
            legend('Turbine Bonete','Spillage Bonete','Turbine Baygorria','Spillage Baygorria',...
                'Turbine Palmar','Turbine Salto Grande','location','southeast');
            xlim([0 max(time)]); ylim([0 1]);
            
            set(0,'CurrentFigure',104); clf(104);
            hold on; grid on; grid minor;
            title('Controls over the Time');
            xlabel('Time');
            ylabel('Control');
            for q = 8:11
                P = plot(time,Controls_OP(:,q));
                P.LineWidth = 1;
            end
            legend(Termicas{1:4});
            xlim([0 max(time)]); ylim([0 1]);

            set(0,'CurrentFigure',105); clf(105); 
            P = plot(time,Controls_OP(:,7)');
            P.LineWidth = 1; hold on;
            title('Controls over the Time');
            xlabel('Time'); ylabel('Control');
            plot(time,zeros(1,length(time)),'r');
            legend('Battery');grid on; grid minor;
            xlim([0 max(time)]); ylim([-0.35 1]);
            
            set(0,'CurrentFigure',106); clf(106);
            plot(time,Vol_OP_1); hold on;
            plot(time,Vol_OP_Real_1);
            ylim([Vmin1 Vmax1]);
            grid on; grid minor;
            xlabel('Normalized Time');
            ylabel('Normalized Volume');
            title('Normalized Volume Bonete');
            legend('Approximated Volume Optimal Path','Real Volume Optimal Path');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',107); clf(107);
            plot(time,Vol_OP_3); hold on;
            plot(time,Vol_OP_Real_3);
            ylim([Vmin3 Vmax3]);
            grid on; grid minor;
            xlabel('Normalized Time');
            ylabel('Normalized Volume');
            title('Normalized Volume Palmar');
            legend('Approximated Volume Optimal Path','Real Volume Optimal Path');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',108); clf(108);
            plot(time,Vol_OP_4); hold on;
            plot(time,Vol_OP_Real_4);
            ylim([Vmin4 Vmax4]);
            grid on; grid minor;
            xlabel('Normalized Time');
            ylabel('Normalized Volume');
            title('Normalized Volume Salto Grande');
            legend('Approximated Volume Optimal Path','Real Volume Optimal Path');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',109); clf(109);
            plot(time,Battery_OP); hold on;
            plot(time,Battery_Real);
            ylim([NAMin NAMax]);
            grid on; grid minor;
            xlabel('Time');
            ylabel('Charge Battery');
            title('Charge Battery');
            legend('Approximated Charge Battery Path','Real Charge Battery Path');
            xlim([0 max(time)]);
                        
            set(0,'CurrentFigure',110); clf(110); hold on;
            P = plot(time,CostControls(1,:)*MWhTokW); P.LineWidth = 1;
            P = plot(time,CostControls(3,:)*MWhTokW); P.LineWidth = 1;
            P = plot(time,CostControls(5,:)*MWhTokW); P.LineWidth = 1;
            P = plot(time,CostControls(6,:)*MWhTokW); P.LineWidth = 1;
            grid on; grid minor;
            title('Instant Controls Costs');
            legend('Bonete Turbined','Baygorria Turbined','Palmar Turbined','Salto Grande Turbined');
            xlabel('Time'); ylabel('USD/MWh');
            
            set(0,'CurrentFigure',111); clf(111); hold on;
            P = plot(time,CostControls(7,:)*MWhTokW);
            P.LineWidth = 1;
            grid on; grid minor;
            title('Instant Controls Costs');
            legend('Battery'); xlabel('Time'); ylabel('USD/MWh');
            
            set(0,'CurrentFigure',112); clf(112);
            hold on;
            for cont = 8:11
                P = plot(time,CostControls(cont,:)*MWhTokW);
                P.LineWidth = 1;
            end
            grid on; grid minor;
            title('Instant Controls Costs');
            legend(Termicas); xlabel('Time'); ylabel('USD/MWh');
            
            set(0,'CurrentFigure',113); clf(113); hold on;
            P = plot(time,CostControls(2,:)); P.LineWidth = 1;
            P = plot(time,CostControls(4,:)); P.LineWidth = 1;
            grid on; grid minor;
            title('Instant Controls Costs');
            legend('Bonete Spillage','Baygorria Spillage');
            xlabel('Time'); ylabel('USD/hm^3');
                        
            set(0,'CurrentFigure',114); clf(114); hold on;
            plot(time,uV1plot/NormCost); plot(time,uV3plot/NormCost);
            plot(time,uV4plot/NormCost); plot(time,uAplot/NormCost);
            grid on; grid minor;
            Leg2 = legend('$u_{V1}$','$u_{V3}$','$u_{V4}$','$u_{A}$','location','southeast');
            set(Leg2,'interpreter', 'latex');
            title('Instant Derivatives Values');
            xlabel('Time');
            ylabel('Derivatives');
            
            set(0,'CurrentFigure',115); clf(115); hold on;
            plot(time,CFL_Plot(:,1));
            grid on; grid minor;
            title('dt/CFL Condition');
            xlabel('Time');
            
            set(0,'CurrentFigure',116); clf(116);
            X = [sum(Power_H1_OP),sum(Power_H2_OP),sum(Power_H3_OP),sum(Power_H4_OP),...
                sum(Power_F1_OP),sum(Power_F2_OP),sum(Power_F3_OP),sum(Power_F4_OP),...
                sum(Power_Battery_OP.*(Power_Battery_OP > 0))];
            txt = {[Dams{1},': '],[Dams{2},': '],[Dams{3},': '],[Dams{4},': '],...
                [Termicas{1},': '],[Termicas{2},': '],[Termicas{3},': '],[Termicas{4},': '],...
                'Battery: '};
            colors = [.5 .5 1;
            .5 .8 1;
            .2 .2 1;
            .8 .8 1;
            1 .6 .6;
            1 .3 .5;
            1 .5 .3;
            1 0 0;
            1 1 .5];
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

            set(0,'CurrentFigure',117); clf(117);
            X = [sum(Cost(1,:)),sum(Cost(2,:)),sum(Cost(3,:)),sum(Cost(4,:)),...
                sum(Cost(5,:)),sum(Cost(6,:)),sum(Cost(7,:)),sum(Cost(8,:))];
            txt = {[Dams{1},': '],[Dams{2},': '],[Dams{3},': '],[Dams{4},': '],...
                [Termicas{1},': '],[Termicas{2},': '],[Termicas{3},': '],[Termicas{4},': ']};
            colors = [.5 .5 1;
            .5 .8 1;
            .2 .2 1;
            .8 .8 1;
            1 .6 .6;
            1 .3 .5;
            1 .5 .3;
            1 0 0];
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
            legend(combinedtxt,'Location','eastoutside');
            colormap(colors);
            
            set(0,'CurrentFigure',118); clf(118); hold on;
            for cont = 7:11
                P = plot(time,CostControls(cont,:)*MWhTokW);
                P.LineWidth = 1;
            end
            grid on; grid minor; legend('Battery',Termicas,'location','southeast')
            title('Instant Controls Costs');
            xlabel('Time'); ylabel('USD/MWh');
            
        end
        
        
        
        toc
    
    end
    
    %% End.

end