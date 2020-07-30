function [HJB_Initial_Cost] = Admissible_Solution_4(OptimalPath,ComputePlots,SavePlots,Expan,Expan_Time,UseLamb21,UseLamb32,Lambda21,Lambda32)

    tic
    
    Phi_V2 = 0.2; % Value of the virtual control for t <= Z21.
    Phi_V3 = 0.2; % Value of the virtual control for t <= Z32.
    
    Derivatives = [1,1,1]; % For debugging.

    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',10000,...
        'FunctionTolerance',1e-10,'MaxIterations',5000,'MaxFunEvals',10^7,'Algorithm','sqp','ConstraintTolerance',1e-7);
    
%     if ComputePlots == 1
%         close all;
%         NeededPlots(1,12,0);
%         if OptimalPath == 1
%             NeededPlots(111,119,0);
%         elseif OptimalPath ~= 0
%             disp('Choose OptimalPath 0 or 1.');
%         return;
%         end
%     elseif ComputePlots ~= 0
%         disp('Choose ComputePlots 0 or 1.');
%         return;
%     end
    
    NT = 2^Expan_Time; % Discretizations in time.
    NV1 = 2^Expan; % Discretizations of Bonete.
    NV2 = 2^Expan; % Discretizations of Baygorria.
    NV3 = 2^Expan; % Discretizations of Palmar.
    NV4 = 2^Expan; % Discretizations of Salto Grande.

    % ==================== Time ====================>>>

    absTime = 0:1/24:1; % Is needed for interpolation.
    Tmax = 1; % Simulation final time.
    Tmin = 0; % Simulation initial time.
    dt = (Tmax-Tmin)/NT;
    time = Tmin:dt:Tmax; % Simulation time.

    % ==================== Parameters (Demand, Water and Fuel) ====================>>

    FlowMax1 = 640; % In m^3/s.
    FlowMax2 = 828; % In m^3/s.
    FlowMax3 = 1373; % In m^3/s.
    FlowMax4 = 4200; % In m^3/s.

    MaxSpill1 = 2420; % In m^3/s.
    MaxSpill2 = 2646; % In m^3/s.
    MaxSpill3 = 2787.5; % In m^3/s.
    MaxSpill4 = 12600; % In m^3/s.

    MaxRel1 = FlowMax1 + MaxSpill1; % In m^3/s.
    MaxRel2 = FlowMax2 + MaxSpill2; % In m^3/s.
    MaxRel3 = FlowMax3 + MaxSpill3; % In m^3/s.
    MaxRel4 = FlowMax4 + MaxSpill4; % In m^3/s.

    MaxV2 = MaxRel1;
    MaxV3 = MaxRel2;

    FuelMax = 2*10^6; % In kW.

    VolMax1 = 9.5*10^9; % In m^3.
    VolMax2 = 617*10^6; % In m^3.
    VolMax3 = 2823*10^6; % In m^3.
    VolMax4 = 5*10^9; % In m^3.

    VolMin1 = 1925*10^6; % In m^3.
    VolMin2 = 424*10^6; % In m^3.
    VolMin3 = 1767*10^6; % In m^3.
    VolMin4 = 3470*10^6; % In m^3.

    TMax = 24*3600; % In s. Number of seconds in a day.

    DMax = 2*10^6; % In kW. The maximum possible value of the demand.

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
    
    if UseLamb21 == 0
        Lambda21 = zeros(1,length(time)+Z21);
    end
    if UseLamb32 == 0
        Lambda32 = zeros(1,length(time)+Z32);
    end
    Lambda21 = interp1(0:1/(length(Lambda21)-1):1,Lambda21,0:1/(length(time)+Z21-1):1);
    Lambda32 = interp1(0:1/(length(Lambda32)-1):1,Lambda32,0:1/(length(time)+Z32-1):1);
    if ComputePlots == 1
        set(0,'CurrentFigure',11);
        plot(0:dt:1+t21,Lambda21); grid on;
        xlim([0 1+t21]);
        set(0,'CurrentFigure',12);
        plot(0:dt:1+t32,Lambda32); grid on;
        xlim([0 1+t32]);
    end

    sigH1 = 0; % Sigma from the SDE of Bonete.
    sigH2 = 0; % Sigma from the SDE of Baygorria.
    sigH3 = 0; % Sigma from the SDE of Palmar.
    sigH4 = 0; % Sigma from the SDE of Salto Grande.

    % Here we use an unreal natural inflow equal to half of the maximum turbine
    % flow, for all time. Also, we compute the maximum and minimum posible
    % volumes of the dams given that natural inputs.
    IT1 = @(t) FlowMax1/3; % In m^3/s.
    IT2 = @(t) FlowMax2/20; % In m^3/s.
    IT3 = @(t) FlowMax3/20; % In m^3/s.
    IT4 = @(t) FlowMax4/3; % In m^3/s.

    Vmax1 = min((Vini1*VolMax1+(FlowMax1+MaxSpill1)*TMax),VolMax1)/VolMax1;
    Vmin1 = max((Vini1*VolMax1-(FlowMax1+MaxSpill1)*TMax),VolMin1)/VolMax1;
    Vmax2 = min((Vini2*VolMax2+(FlowMax2+MaxSpill2)*TMax),VolMax2)/VolMax2;
    Vmin2 = max((Vini2*VolMax2-(FlowMax2+MaxSpill2)*TMax),VolMin2)/VolMax2;
    Vmax3 = min((Vini3*VolMax3+(FlowMax3)*TMax),VolMax3)/VolMax3;
    Vmin3 = max((Vini3*VolMax3-(FlowMax3)*TMax),VolMin3)/VolMax3;
    Vmax4 = min((Vini4*VolMax4+(FlowMax4)*TMax),VolMax4)/VolMax4;
    Vmin4 = max((Vini4*VolMax4-(FlowMax4)*TMax),VolMin4)/VolMax4;
    
    RedMax1 = Vmax1-Vmin1;
    RedMax2 = Vmax2-Vmin2;
    RedMax3 = Vmax3-Vmin3;
    RedMax4 = Vmax4-Vmin4;
    Vmax1 = Vmax1 - RedMax1/3;
    Vmin1 = Vmin1 + RedMax1/3;
    Vmax2 = Vmax2 - RedMax2/4;
    Vmin2 = Vmin2 + RedMax2/4;
    Vmax3 = Vmax3 - RedMax3/3;
    Vmin3 = Vmin3 + RedMax3*0;
    Vmax4 = Vmax4 - RedMax4/3;
    Vmin4 = Vmin4 + RedMax4*0;
    
    load('D01012017.mat'); % Demand from 01/01/2017. It loads a 25x1 vector named D01012017.

    DV1 = Vmax1-Vmin1;
    DV2 = Vmax2-Vmin2;
    DV3 = Vmax3-Vmin3;
    DV4 = Vmax4-Vmin4;

    dV1 = DV1/NV1;
    dV2 = DV2/NV2;
    dV3 = DV3/NV3;
    dV4 = DV4/NV4;

    V1 = Vmin1:dV1:Vmax1; % Discretized volume of Bonete.
    V2 = Vmin2:dV2:Vmax2; % Discretized volume of Baygorria.
    V3 = Vmin3:dV3:Vmax3; % Discretized volume of Palmar.
    V4 = Vmin4:dV4:Vmax4; % Discretized volume of Salto Grande.

    [vplot1_4,vplot3_4] = meshgrid(V1,V3); % For plotting surfaces.
    [vplot1_3,vplot4_3] = meshgrid(V1,V4); % For plotting surfaces.
    [vplot3_1,vplot4_1] = meshgrid(V3,V4); % For plotting surfaces.

%     D = DMax*(0.5-(0.4)*sin(2*pi*time)); % In kW. (THIS (1))

%     D = DMax*(0.5)*time./time; % In kW. (THIS (2))

    D = D01012017'*1000*1.2; % In kW. (THIS (3))
    D = interp1(absTime,D,time)+DMax*0.15*sin(4*pi*time); % Interpolation. (THIS (3))

    % ==================== Plot Demand (plot 6) ====================>>

    if ComputePlots == 1
        set(0,'CurrentFigure',1);
        hold on;
        xlabel('Normalized Time');
        title('Normalized Demand');
        grid on;
        plot(time,D/DMax);
        ylim([0 1]);
        pause(0.001);
    end

    % ==================== Matrices ====================>>

    u = zeros(length(V1),length(V3),length(V4));
    HAM = zeros(length(V1),length(V3),length(V4));
    U = zeros(length(V1),length(V3),length(V4),length(time));
    uV1 = zeros(length(V1),length(V3),length(V4));
    uV3 = zeros(length(V1),length(V3),length(V4));
    uV4 = zeros(length(V1),length(V3),length(V4));
    UV = {};
    ControlsFMC = zeros(length(V1),length(V3),length(V4),9); % F,T1,S1,T2,S2,T3,T4,V2,V3.
    ControlsFMC_Aux = ControlsFMC;
    Controls_T = ControlsFMC;

    % ==================== More parameters (Water and Fuel) ====================>>

    H1 = @(V1) (-3.74)*(V1).^2 + (16.7)*(V1) + 67.7; % Height of the water as a function of the normal volume (Bonete).
    % H2 = @(V2) (-3.25)*(V2).^2 + (12.1)*(V2) + 45.6; % Height of the water as a function of the normal volume (Baygorria).
    H2 = 54; % We fix the level of Baygorria to its nominal level.
    H3 = @(V3) (-3.76)*(V3).^2 + (16.9)*(V3) + 26.8; % Height of the water as a function of the normal volume (Palmar).
    H4 = @(V4) (-19.8)*(V4).^2 + (51.5)*(V4) + 3.79; % Height of the water (m) as a function of the normal volume (Salto Grande).

    cT1 = 1.02;
    cS1 = cT1*MaxSpill1/FlowMax1;
    d1 = @(tur,spill) cT1*tur+cS1*spill; % Variation in the high due of the released flow (Bonete).

    cT2 = 0.57;
    cS2 = cT2*MaxSpill2/FlowMax2;
    d2 = @(tur,spill) cT2*tur+cS2*spill; % Variation in the high due of the released flow (Baygorria).

    cT3 = 3.77;
    cS3 = cT3*MaxSpill3/FlowMax3;
    d3 = @(tur,spill) cT3*tur+cS3*spill; % Variation in the high due of the released flow (Palmar).

    cT4 = 5.34;
    cS4 = cT4*MaxSpill4/FlowMax4;
    d4 = @(tur,spill) cT4*tur+cS4*spill; % Variation in the high due of the released flow (Salto Grande).

    eta = 8.72; % Efficiency.

    h04 = 5.1; % Level after Salto Grande. Now is a constant, but may depend of the time.
    h03 = 5.1; % Level after Palmar. Now is a constant, but may depend of the time.

    Khe1 = 0.3; % In Mil USD/MWh. Cost per energy of Bonete.
    Khe2 = 0; % In Mil USD/MWh.
    Khe3 = 0.3; % In Mil USD/MWh.
    Khe4 = 0.3; % In Mil USD/MWh.
    Kf = 13; % In Mil USD/MWh. Cost per energy of the fuel.

    CH1 = Khe1*eta*(H1(Vini1)-d1(1,0)-H2)/(3.6*10^6); % Water's value of Bonete.
    CH2 = Khe2*eta*(H2-d2(1,0)-H3(Vini3))/(3.6*10^6); % Water's value of Baygorria.
    CH3 = Khe3*eta*(H3(Vini3)-d3(1,0)-h03)/(3.6*10^6); % Water's value of Palmar.
    CH4 = Khe4*eta*(H4(Vini4)-d4(1,0)-h04)/(3.6*10^6); % Water's value of Salto Grande.
    CF = Kf/(3.6*10^6);

    % We divide over 3600 to pass from hours to seconds and over 1000 to move
    % from Mega to Kilo.

    olCF = CF*FuelMax;
    olCT1 = @(uV1,t) FlowMax1*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCS1 = @(uV1,t) MaxSpill1*(CH1 - Lambda21(t+Z21) - uV1/VolMax1);
    olCT2 = @(uV2,t) FlowMax2*(CH2 - Lambda32(t+Z32) - uV2/VolMax2 + Lambda21(t));
    olCS2 = @(uV2,t) MaxSpill2*(CH2 - Lambda32(t+Z32) - uV2/VolMax2 + Lambda21(t));
    olCT3 = @(uV3,t) FlowMax3*(CH3 - uV3/VolMax3);
    olCT4 = @(uV4,t) FlowMax4*(CH4 - uV4/VolMax4);
    H0 = @(uV1,uV2,uV3,uV4,t) IT1(t)*uV1/VolMax1 + IT3(t)*uV3/VolMax3 + IT4(t)*uV4/VolMax4 - Lambda21(t)*IT2(t);
    olCV3 = @(uV3,t) MaxV3*(Lambda32(t) + uV3/VolMax3);

    d = @(uV1,uV2,uV3,uV4,t) [olCF,olCT1(uV1,t),olCS1(uV1,t),olCT2(uV2,t),olCS2(uV2,t),olCT3(uV3,t),olCT4(uV4,t)]';

    k1 = @(t) IT2(t)/MaxV2;
    k2 = FlowMax2/MaxV2;
    k3 = MaxSpill2/MaxV2;
    k = @(t) [k1(t),k2,k3];

    % Hydraulic power: PH = T^2*K2 + T*K1 + T*S*K3.

    K11 = @(V1) eta*FlowMax1*(H1(V1)-H2);
    K12 = @(V3) eta*FlowMax2*(H2-H3(V3));
    K13 = @(V3) eta*FlowMax3*(H3(V3)-h03);
    K14 = @(V4) eta*FlowMax4*(H4(V4)-h04);

    K21 = -eta*FlowMax1*cT1;
    K22 = -eta*FlowMax2*cT2;
    K23 = -eta*FlowMax3*cT3;
    K24 = -eta*FlowMax4*cT4;

    K31 = -eta*FlowMax1*cS1;
    K32 = -eta*FlowMax2*cS2;
    K33 = -eta*FlowMax3*cS3;
    K34 = -eta*FlowMax4*cS4;

    b = @(V1,V3,V4) [FuelMax,K11(V1),0,K12(V3),0,K13(V3),K14(V4)]';

    Q = [0 0 0 0 0 0 0;
        0 K21 K31/2 0 0 0 0;
        0 K31/2 0 0 0 0 0;
        0 0 0 K22 K32/2 0 0;
        0 0 0 K32/2 0 0 0;
        0 0 0 0 0 K23 0;
        0 0 0 0 0 0 K24]; % Is fixed.

    c = @(t) -D(t);

    % ==================== Simulation ====================>>

    for t = length(time):-1:1

        % ==================== Finite Differences Section ====================>>

        for v1 = 1:length(V1)

            for v3 = 1:length(V3)

                for v4 = 1:length(V4)

                    if (v1~=1) && (v3~=1) && (v4~=1) && (v1~=length(V1)) && (v3~=length(V3)) && (v4~=length(V4)) && (t<length(time)-1)

                        if (IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0) && Derivatives(1)
                            uV1(v1,v3,v4) = (u(v1+1,v3,v4)-u(v1,v3,v4))/(dV1);
                        elseif Derivatives(1)
                            uV1(v1,v3,v4) = (u(v1,v3,v4)-u(v1-1,v3,v4))/(dV1);
                        end

                        if (IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,9) >= 0) && Derivatives(2)
                            uV3(v1,v3,v4) = (u(v1,v3+1,v4)-u(v1,v3,v4))/(dV3);
                        elseif Derivatives(2)
                            uV3(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3-1,v4))/(dV3);
                        end

                        if (IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0) && Derivatives(3)
                            uV4(v1,v3,v4) = (u(v1,v3,v4+1)-u(v1,v3,v4))/(dV4);
                        elseif Derivatives(3)
                            uV4(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3,v4-1))/(dV4);
                        end

                    elseif t < length(time) - 1

                        if v1 == 1 && Derivatives(1)
                            if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                                uV1(v1,v3,v4) = (u(v1+2,v3,v4)-u(v1+1,v3,v4))/(dV1);
                            else
                                uV1(v1,v3,v4) = (u(v1+1,v3,v4)-u(v1,v3,v4))/(dV1);
                            end
                        elseif v1 == length(V1) && Derivatives(1)
                            if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                                uV1(v1,v3,v4) = (u(v1,v3,v4)-u(v1-1,v3,v4))/(dV1);
                            else
                                uV1(v1,v3,v4) = (u(v1-1,v3,v4)-u(v1-2,v3,v4))/(dV1);
                            end
                        elseif (v3==1 || v3==length(V3) || v4==1 || v4==length(V4)) && (v1~=1) && (v1~=length(V1)) && Derivatives(1)
                            if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                                uV1(v1,v3,v4) = (u(v1+1,v3,v4)-u(v1,v3,v4))/(dV1);
                            else
                                uV1(v1,v3,v4) = (u(v1,v3,v4)-u(v1-1,v3,v4))/(dV1);
                            end
                        end

                        if v3 == 1 && Derivatives(2)
                            if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,9) >= 0
                                uV3(v1,v3,v4) = (u(v1,v3+2,v4)-u(v1,v3+1,v4))/(dV3);
                            else
                                uV3(v1,v3,v4) = (u(v1,v3+1,v4)-u(v1,v3,v4))/(dV3);
                            end
                        elseif v3 == length(V3) && Derivatives(2)
                            if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,9) >= 0
                                uV3(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3-1,v4))/(dV3);
                            else
                                uV3(v1,v3,v4) = (u(v1,v3-1,v4)-u(v1,v3-2,v4))/(dV3);
                            end
                        elseif (v1==1 || v1==length(V1) || v4==1 || v4==length(V4)) && (v3~=1) && (v3~=length(V3)) && Derivatives(2)
                            if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,9) >= 0
                                uV3(v1,v3,v4) = (u(v1,v3+1,v4)-u(v1,v3,v4))/(dV3);
                            else
                                uV3(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3-1,v4))/(dV3);
                            end
                        end

                        if v4 == 1 && Derivatives(3)
                            if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0
                                uV4(v1,v3,v4) = (u(v1,v3,v4+2)-u(v1,v3,v4+1))/(dV4);
                            else
                                uV4(v1,v3,v4) = (u(v1,v3,v4+1)-u(v1,v3,v4))/(dV4);
                            end
                        elseif v4 == length(V4) && Derivatives(3)
                            if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0
                                uV4(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3,v4-1))/(dV4);
                            else
                                uV4(v1,v3,v4) = (u(v1,v3,v4-1)-u(v1,v3,v4-2))/(dV4);
                            end
                        elseif (v1==1 || v1==length(V1) || v3==1 || v3==length(V3)) && (v4~=1) && (v4~=length(V4)) && Derivatives(3)
                            if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0
                                uV4(v1,v3,v4) = (u(v1,v3,v4+1)-u(v1,v3,v4))/(dV4);
                            else
                                uV4(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3,v4-1))/(dV4);
                            end
                        end

                    end

                end

            end

        end
        
        UV{t} = {uV1,uV3,uV4};

        % ==================== Controls Section ====================>>

        ControlsFMC_Aux = Controls_T;
        
        TL = (length(V1)*length(V3)*length(V4));
        LV1 = length(V1);
        LV3 = length(V3);
        LV4 = length(V4);
        Size_Mat = size(u);
        uV1_PF = uV1(:);
        uV3_PF = uV3(:);
        uV4_PF = uV4(:);
        Controls_Aux = zeros(TL,9);
        Warm_Start = reshape(Controls_T,LV1,LV3,LV4,9);
        Warm_Start = Warm_Start(:,:,:,1:7);
        
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
            else % if UseLamb21 == 0
                Aux_V2 = 0;
                PhiV2Fixed = 1;
            end
            
            [ind_V1,ind_V3,ind_V4] = ind2sub(Size_Mat,index);
            WS = squeeze(Warm_Start(ind_V1,ind_V3,ind_V4,:));
%             WS = zeros(7,1);
        
            Controls_Aux(index,:) = [Quad_FMC_Lambda(Q,b(V1(ind_V1),V3(ind_V3),V4(ind_V4)),c(t),...
            	d(uV1_PF(index),0,uV3_PF(index),uV4_PF(index),t),...
            	k(t),options,PhiV2Fixed,Aux_V2,WS);Aux_V3];
        
        end
        
        Controls_T = reshape(Controls_Aux,LV1,LV3,LV4,9);
                
        % ==================== Time Step ====================>>

        if t ~= length(time)

            for v1 = 1:length(V1)

                for v3 = 1:length(V3)

                    for v4 = 1:length(V4)

                        Controls(1:9) = ControlsFMC_Aux(v1,v3,v4,:);
                        
                        HAM(v1,v3,v4) = TMax*(Controls*...
                            [d(uV1(v1,v3,v4),0,uV3(v1,v3,v4),uV4(v1,v3,v4),t);0;olCV3(uV3(v1,v3,v4),t)]+...
                            H0(uV1(v1,v3,v4),0,uV3(v1,v3,v4),uV4(v1,v3,v4),t));

                        u(v1,v3,v4) = u(v1,v3,v4) + dt*HAM(v1,v3,v4);

                    end

                end

            end
            
            U(:,:,:,t) = u;

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
            clf(2);
            hold on;
            for index = 1:floor(length(V4)/2):length(V4)
                S = surf(vplot1_4,vplot3_4,squeeze(u(:,:,index)/1e5)');
                S.EdgeColor = 'none';
            end
            xlabel('Vol. Bonete');
            ylabel('Vol. Palmar');
            zlabel('Normalized Cost');
            title(['Optimal Cost Function at Time = ',num2str(t)]);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin3 Vmax3]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V4_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',3);
            clf(3);
            hold on;
            for index = 1:floor(length(V4)/2):length(V4)
                S = surf(vplot1_3,vplot4_3,squeeze(u(:,index,:)/1e5)');
                S.EdgeColor = 'none';
            end
            xlabel('Vol. Bonete');
            ylabel('Vol. Salto Grande');
            zlabel('Normalized Cost');
            title(['Optimal Cost Function at Time = ',num2str(t)]);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin4 Vmax4]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V3_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',4);
            clf(4);
            hold on;
            for index = 1:floor(length(V4)/2):length(V4)
                S = surf(vplot3_1,vplot4_1,squeeze(u(index,:,:)/1e5)');
                S.EdgeColor = 'none';
            end
            xlabel('Vol. Palmar');
            ylabel('Vol. Salto Grande');
            zlabel('Normalized Cost');
            title(['Optimal Cost Function at Time = ',num2str(t)]);
            xlim([Vmin3 Vmax3]);
            ylim([Vmin4 Vmax4]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V1_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',5);
            clf(5);
            hold on;
%             for index=1:1:length(V4)
                S = surf(vplot1_4,vplot3_4,squeeze(HAM(:,:,index))');
                S.EdgeColor = 'none';
%             end
            xlabel('Vol. Bonete');
            ylabel('Vol. Palmar');
            title(['Hamiltonian at time = ',num2str(t)]);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin3 Vmax3]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V4_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',6);
            clf(6);
            hold on;
%             for index=1:1:length(V3)
                S = surf(vplot1_3,vplot4_3,squeeze(HAM(:,index,:))');
                S.EdgeColor = 'none';
%             end
            xlabel('Vol. Bonete');
            ylabel('Vol. Salto Grande');
            title(['Hamiltonian at time = ',num2str(t)]);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin4 Vmax4]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V3_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',7);
            clf(7);
            hold on;
%             for index=1:1:length(V1)
                S = surf(vplot3_1,vplot4_1,squeeze(HAM(index,:,:))');
                S.EdgeColor = 'none';
%             end
            xlabel('Vol. Palmar');
            ylabel('Vol. Salto Grande');
            title(['Hamiltonian at time = ',num2str(t)]);
            xlim([Vmin3 Vmax3]);
            ylim([Vmin4 Vmax4]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V1_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',8);
            clf(8);
            hold on;
%             for index=1:1:length(V4)
                S = surf(vplot1_4,vplot3_4,squeeze(uV4(:,:,index)/1e5)');
                S.EdgeColor = 'none';
%             end
            xlabel('Vol. Bonete');
            ylabel('Vol. Palmar');
            title(['u_{V4} at Time = ',num2str(t)]);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin3 Vmax3]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V4_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',9);
            clf(9);
            hold on;
%             for index=1:1:length(V3)
                S = surf(vplot1_3,vplot4_3,squeeze(uV3(:,index,:)/1e5)');
                S.EdgeColor = 'none';
%             end
            xlabel('Vol. Bonete');
            ylabel('Vol. Salto Grande');
            title(['u_{V3} at Time = ',num2str(t)]);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin4 Vmax4]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V3_',num2str(Expan)],'epsc');
            end
            
            set(0,'CurrentFigure',10);
            clf(10);
            hold on;
%             for index=1:1:length(V1)
                S = surf(vplot3_1,vplot4_1,squeeze(uV1(index,:,:)/1e5)');
                S.EdgeColor = 'none';
%             end
            xlabel('Vol. Palmar');
            ylabel('Vol. Salto Grande');
            title(['u_{V1} at Time = ',num2str(t)]);
            xlim([Vmin3 Vmax3]);
            ylim([Vmin4 Vmax4]);
            zlim auto; box; view(-40,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_V1_',num2str(Expan)],'epsc');
            end

        elseif ComputePlots ~= 0
            disp('Choose ComputePlots between 0 or 1.');
            return;
        end    
        
        fprintf('Discretizations = %d. Completed : %d',NV4,100*(1-t/length(time)));
        fprintf('\n');
        pause(0.1);

    end
    
    IP1 = min(abs(V1-Vini1)) == abs(V1-Vini1);
    IP3 = min(abs(V3-Vini3)) == abs(V3-Vini3);
    IP4 = min(abs(V4-Vini4)) == abs(V4-Vini4);
    HJB_Initial_Cost = u(IP1,IP3,IP4);

    toc
    
    %% Optimal Path Section:
    
    if OptimalPath == 1

        tic
        
        IP1 = find(min(abs(V1-Vini1)) == abs(V1-Vini1));
        IP3 = find(min(abs(V3-Vini3)) == abs(V3-Vini3));
        IP4 = find(min(abs(V4-Vini4)) == abs(V4-Vini4));
        HJB_Initial_Cost = u(IP1,IP3,IP4);

        Optimal_Cost = zeros(1,length(time));
        Phi_V2_Aux = zeros(1,length(time)); % To save the past virtual control V2.
        Phi_V3_Aux = zeros(1,length(time)); % To save the past virtual control V3.
        Controls_OP = zeros(length(time),9);
        Power_F_OP = zeros(1,length(time));
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
        Cost = zeros(5,length(time));
        CostControls = zeros(7,length(time));

        x0 = zeros(7,1);
        x0(1:7) = ControlsFMC(IP1,IP3,IP4,1:7);

        Power_F = @(x1) x1*FuelMax; % In kW.
        Power_H1 = @(x2,x3,V1) eta*x2*FlowMax1*(H1(V1)-d1(x2,x3)-H2); % In kW.
        Power_H2 = @(x4,x5,V3) eta*x4*FlowMax2*(H2-d2(x4,x5)-H3(V3)); % In kW.
        Power_H3 = @(x6,V3) eta*x6*FlowMax3*(H3(V3)-d3(x6,0)-h03); % In kW.
        Power_H4 = @(x7,V4) eta*x7*FlowMax4*(H4(V4)-d4(x7,0)-h04); % In kW.
        
        red_dt = 1;
        newdt = dt*red_dt;
        rest = length(time)-1;
        while rest > 0
            rest = rest - red_dt;
        end
        Z21 = Z21-mod(Z21,red_dt);
        Z32 = Z32-mod(Z32,red_dt);
        
        for t = 1:red_dt:length(time) + rest
            
            uV1 = UV{t}{1}(IP1,IP3,IP4); uV3 = UV{t}{2}(IP1,IP3,IP4); uV4 = UV{t}{3}(IP1,IP3,IP4); uV2 = 0;

            % ==================== Controls Section ====================>>

            if t > Z21
                if UseLamb21 == 1
                    Controls_OP(t,1:8) = Quad_FMC_Lambda(Q,b(Vol_OP_1(t),Vol_OP_3(t),Vol_OP_4(t)),c(t),...
                        d(uV1,0,uV3,uV4,t),k(t),options,1,Phi_V2_Aux(t-Z21),x0);
                else
                    Controls_OP(t,1:8) = Quad_FMC_Lambda(Q,b(Vol_OP_1(t),Vol_OP_3(t),Vol_OP_4(t)),c(t),...
                        d(uV1,0,uV3,uV4,t),k(t),options,1,0,x0);
                end
            else
                if UseLamb21 == 1
                    Controls_OP(t,1:8) = Quad_FMC_Lambda(Q,b(Vol_OP_1(t),Vol_OP_3(t),Vol_OP_4(t)),c(t),...
                        d(uV1,0,uV3,uV4,t),k(t),options,1,Phi_V2,x0);
                else
                    Controls_OP(t,1:8) = Quad_FMC_Lambda(Q,b(Vol_OP_1(t),Vol_OP_3(t),Vol_OP_4(t)),c(t),...
                        d(uV1,0,uV3,uV4,t),k(t),options,1,0,x0);
                end
            end

            Phi_V2_Aux(t) = (Controls_OP(t,2)*FlowMax1 + Controls_OP(t,3)*MaxSpill1) / MaxV2;
            Phi_V3_Aux(t) = (Controls_OP(t,4)*FlowMax2 + Controls_OP(t,5)*MaxSpill2) / MaxV3;
            x0(1:7) = Controls_OP(t,1:7);
            
            if UseLamb32 == 1
                if t > Z32
                    Controls_OP(t,9) = Phi_V3_Aux(t-Z32);
                else
                    Controls_OP(t,9) = Phi_V3;
                end
            else
                Controls_OP(t,9) = 0;
            end
            
            Aux_Cost = Controls_OP(t,1:7).*d(uV1,0,uV3,uV4,t);
            
            Cost(1,t) = Aux_Cost(1);
            Cost(2,t) = Aux_Cost(2)+Aux_Cost(3);
            Cost(3,t) = Aux_Cost(4)+Aux_Cost(5);
            Cost(4,t) = Aux_Cost(6);
            Cost(5,t) = Aux_Cost(7);
            
            CostControls(:,t) = d(uV1,0,uV3,uV4,t);

            % ==================== Dynamics and Powers Section ====================>>

            if t ~= length(time) + rest
                
                if t+red_dt <= length(time)
                    
                    Vol_OP_Real_1(t+red_dt) = Vol_OP_Real_1(t) + (newdt*TMax/VolMax1)*(IT1(t)-...
                        Controls_OP(t,2)*FlowMax1-Controls_OP(t,3)*MaxSpill1);
                    Vol_OP_Real_3(t+red_dt) = Vol_OP_Real_3(t) + (newdt*TMax/VolMax3)*(IT3(t)-...
                        Controls_OP(t,6)*FlowMax3+Controls_OP(t,9)*MaxV3);
                    Vol_OP_Real_4(t+red_dt) = Vol_OP_Real_4(t) + (newdt*TMax/VolMax4)*(IT4(t)-...
                        Controls_OP(t,7)*FlowMax4);

                    % This is to approximate with the point of the grid.

                    AuxV1 = V1(min(abs(V1-Vol_OP_Real_1(t+red_dt))) == abs(V1-Vol_OP_Real_1(t+red_dt)));
                    Vol_OP_1(t+red_dt) = AuxV1(1);
                    AuxV3 = V3(min(abs(V3-Vol_OP_Real_3(t+red_dt))) == abs(V3-Vol_OP_Real_3(t+red_dt)));
                    Vol_OP_3(t+red_dt) = AuxV3(1);
                    AuxV4 = V4(min(abs(V4-Vol_OP_Real_4(t+red_dt))) == abs(V4-Vol_OP_Real_4(t+red_dt)));
                    Vol_OP_4(t+red_dt) = AuxV4(1);

                    IP1 = find(V1 == Vol_OP_1(t+red_dt));
                    IP3 = find(V3 == Vol_OP_3(t+red_dt));
                    IP4 = find(V4 == Vol_OP_4(t+red_dt));

                    Optimal_Cost(t+red_dt) = Optimal_Cost(t) + newdt*TMax*(Controls_OP(t,:)*[d(uV1,0,uV3,uV4,t);0;olCV3(uV3,t)] + H0(uV1,uV2,uV3,uV4,t));
                
                end
                
            end

            % ==================== 'Others' Section ====================>>

            Power_F_OP(t) = Power_F(Controls_OP(t,1));
            Power_H1_OP(t) = Power_H1(Controls_OP(t,2),Controls_OP(t,3),Vol_OP_1(t));
            Power_H2_OP(t) = Power_H2(Controls_OP(t,4),Controls_OP(t,5),Vol_OP_3(t));
            Power_H3_OP(t) = Power_H3(Controls_OP(t,6),Vol_OP_3(t));
            Power_H4_OP(t) = Power_H4(Controls_OP(t,7),Vol_OP_4(t));

        end

        if ComputePlots == 1
            
            if red_dt > 1
                for j = length(time):-1:1

                    if Vol_OP_Real_1(j) == 0

                        Vol_OP_1(j) = [];
                        Vol_OP_3(j) = [];
                        Vol_OP_4(j) = [];
                        Vol_OP_Real_1(j) = [];
                        Vol_OP_Real_3(j) = [];
                        Vol_OP_Real_4(j) = [];
                        time(j) = [];
                        Power_H2_OP(j) = [];
                        Power_H1_OP(j) = [];
                        Power_H3_OP(j) = [];
                        Power_H4_OP(j) = [];
                        Power_F_OP(j) = [];
                        D(j) = [];
                        Optimal_Cost(j) = [];
                        Controls_OP(j,:) = [];

                    end

                end
            end

            set(0,'CurrentFigure',111);
            grid on;
            hold on;
            area(time,[Power_H2_OP;Power_H1_OP;Power_H3_OP;Power_H4_OP;Power_F_OP]')
            [~,h_legend] = legend('Baygorria','Bonete','Palmar','Salto Grande','Fuel');
            plot(time,D,'*');
            xlabel('Normalized Time');
            title('Power Balance')
            ylabel('Power (kW)');
            xlim([0 max(time)]);

            set(0,'CurrentFigure',112);
            hold on;
            plot(time,Optimal_Cost/1e5);
            plot(time,ones(1,length(time))*HJB_Initial_Cost/1e5);
            grid on;
            legend('Accumulated Cost','HJB Final Cost','location','southeast');
            xlabel('Normalized Time');
            ylabel('USD');
            title('Costs Comparison');
            xlim([0 max(time)]);

            set(0,'CurrentFigure',113);
            ylim([0,1]);
            hold on;
            grid on;
            title('Controls over the Time');
            xlabel('Normalized Time');
            for q = 1:7
                plot(time,Controls_OP(:,q));
            end
            legend('Fuel','Turbine Bonete','Spillage Bonete','Turbine Baygorria','Spillage Baygorria',...
                'Turbine Palmar','Turbine Salto Grande');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',114);
            plot(time,Vol_OP_1);
            hold on;
            plot(time,Vol_OP_Real_1);
            ylim([Vmin1 Vmax1]);
            grid on;
            xlabel('Normalized Time');
            title('Normalized Volume Bonete');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',115);
            plot(time,Vol_OP_3);
            hold on;
            plot(time,Vol_OP_Real_3);
            ylim([Vmin3 Vmax3]);
            grid on;
            xlabel('Normalized Time');
            title('Normalized Volume Palmar');
            xlim([0 max(time)]);
            
            set(0,'CurrentFigure',116);
            plot(time,Vol_OP_4);
            hold on;
            plot(time,Vol_OP_Real_4);
            ylim([Vmin4 Vmax4]);
            grid on;
            xlabel('Normalized Time');
            title('Normalized Volume Salto Grande');
            xlim([0 max(time)]);
            
            if UseLamb21 == 0 && UseLamb32 == 0
            
                set(0,'CurrentFigure',117);
                X = [sum(Power_H1_OP),sum(Power_H3_OP),sum(Power_H4_OP),sum(Power_F_OP)];
                p = pie(X);
                title('Energy Distribution');
                pText = findobj(p,'Type','text');
                percentValues = get(pText,'String');
                txt = {'Bonete: ';'Palmar: ';'Salto Grande: ';'Fuel: '}; 
                combinedtxt = strcat(txt,percentValues);
                pText(1).String = combinedtxt(1);
                pText(2).String = combinedtxt(2);
                pText(3).String = combinedtxt(3);
                pText(4).String = combinedtxt(4);
                colormap([.5 .8 1; % Bonete.
                .5 .5 1; % Palmar.
                0 0 1; % Salto Grande.
                1 .5 .5]); % Fuel.

                set(0,'CurrentFigure',118);
                X = [sum(Cost(2,:)),sum(Cost(4,:)),sum(Cost(5,:)),sum(Cost(1,:))];
                p = pie(X);
                title('Costs Distribution');
                pText = findobj(p,'Type','text');
                percentValues = get(pText,'String'); 
                txt = {' ';' ';'All Dams Cost: ';'Fuel Cost: '}; 
                combinedtxt = strcat(txt,percentValues);
                pText(1).String = ' ';
                pText(2).String = ' ';
                pText(3).String = combinedtxt(3);
                pText(4).String = combinedtxt(4);
                colormap([.5 .8 1; % Bonete.
                .5 .5 1; % Palmar.
                0 0 1; % Salto Grande.
                1 .5 .5]); % Fuel.

            elseif UseLamb21 == 1 && UseLamb32 == 1
                
                set(0,'CurrentFigure',117);
                X = [sum(Power_H1_OP),sum(Power_H2_OP),sum(Power_H3_OP),sum(Power_H4_OP),sum(Power_F_OP)];
                p = pie(X);
                title('Energy Distribution');
                pText = findobj(p,'Type','text');
                percentValues = get(pText,'String');
                txt = {'Bonete: ';'Baygorria: ';'Palmar: ';'Salto Grande: ';'Fuel: '}; 
                combinedtxt = strcat(txt,percentValues);
                pText(1).String = combinedtxt(1);
                pText(2).String = combinedtxt(2);
                pText(3).String = combinedtxt(3);
                pText(4).String = combinedtxt(4);
                pText(5).String = combinedtxt(5);
                colormap([.5 .8 1; % Bonete.
                .5 .5 1; % Baygorria.
                .3 .3 1; % Palmar.
                0 0 1; % Salto Grande.
                1 .5 .5]); % Fuel.

                set(0,'CurrentFigure',118);
                X = [sum(Cost(2,:)),sum(Cost(3,:)),sum(Cost(4,:)),sum(Cost(5,:)),sum(Cost(1,:))];
                p = pie(X);
                title('Costs Distribution');
                pText = findobj(p,'Type','text');
                percentValues = get(pText,'String'); 
                if sum(Cost(3,:)) <= 0
                    txt = {' ';' ';'All Dams Cost: ';'Fuel Cost: '};
                    combinedtxt = strcat(txt,percentValues);
                    pText(1).String = ' ';
                    pText(2).String = ' ';
                    pText(3).String = combinedtxt(3);
                    pText(4).String = combinedtxt(4);
                    colormap([.5 .8 1; % Bonete.
                    .3 .3 1; % Palmar.
                    0 0 1; % Salto Grande.
                    1 .5 .5]); % Fuel.
                else
                    txt = {' ';' ';'All Dams Cost: ';' ';'Fuel Cost: '};
                    combinedtxt = strcat(txt,percentValues);
                    pText(1).String = ' ';
                    pText(2).String = ' ';
                    pText(3).String = combinedtxt(3);
                    pText(4).String = ' ';
                    pText(5).String = combinedtxt(5);
                    colormap([.5 .8 1; % Bonete.
                    .5 .5 1; % Baygorria.
                    .3 .3 1; % Palmar.
                    0 0 1; % Salto Grande.
                    1 .5 .5]); % Fuel.
                end

            end
            
            set(0,'CurrentFigure',119);
            hold on;
            for cont = 2:7
                plot(time,CostControls(cont,:));
                
            end
            grid on;
            legend('Bonete Turbined','Bonete Spillage','Baygorria Turbined','Baygorria Spillage',...
                'Palmar Turbined','Salto Grande Turbined');
            xlabel('Normalized Time');
            ylabel('Normalized Price');
            
        end

        toc
    
    end

end