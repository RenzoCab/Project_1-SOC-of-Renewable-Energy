function [Uini] = Cost_Ini_Deterministic(ComputePlots,SavePlots,expansion)

    close all

    Lambda21 = zeros(1,1000);
    Lambda32 = zeros(1,1000);
    options = optimoptions('fmincon','Display','notify','MaxFunctionEvaluations',10000,...
        'FunctionTolerance',1.0000e-5,'MaxIterations',3000,'MaxFunEvals',10^6,'Algorithm','sqp');

    NeededPlots(2);
    NT = 2*2^expansion; % Discretizations in time.
    NV1 = 2*2^expansion; % Discretizations of Bonete.
    NV2 = 2*2^expansion; % Discretizations of Baygorria.
    NV3 = 2*2^expansion; % Discretizations of Palmar.
    NV4 = 2*2^expansion; % Discretizations of Salto Grande.

    % ==================== Time ====================>>>

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
    Vini3 = 0.9; % Initial condition between 0.5 and 1.
    Vini4 = 0.9; % Initial condition between 0.6 and 1.

    t21 = 1/4; % Delay between Bonete and Baygorria.
    Z21 = floor(length(time)*t21); % Discrete delay between Bonete and Baygorria.
    t32 = 1/6; % Delay between Baygorria and Palmar.
    Z32 = floor(length(time)*t32); % Discrete delay between Baygorria and Palmar.
    t31 = 10/24; % Delay between Bonete and Palmar.
    Z31 = floor(length(time)*t31); % Discrete delay between Bonete and Palmar.

    sigH1 = 0; % Sigma from the SDE of Bonete.
    sigH2 = 0; % Sigma from the SDE of Baygorria.
    sigH3 = 0; % Sigma from the SDE of Palmar.
    sigH4 = 0; % Sigma from the SDE of Salto Grande.

    % Here we use an unreal natural inflow equal to half of the maximum turbine
    % flow, for all time. Also, we compute the maximum and minimum posible
    % volumes of the dams given that natural inputs.
    IT1 = @(t) FlowMax1/2; % In m^3/s.
    IT2 = @(t) FlowMax2/2; % In m^3/s.
    IT3 = @(t) FlowMax3/2; % In m^3/s.
    IT4 = @(t) FlowMax4/2; % In m^3/s.
    Vmax1 = min((Vini1*VolMax1+IT1(1)*TMax),VolMax1)/VolMax1;
    Vmin1 = max((Vini1*VolMax1-(FlowMax1+MaxSpill1)*TMax),VolMin1)/VolMax1;
    Vmax2 = min((Vini2*VolMax2+IT2(1)*TMax),VolMax2)/VolMax2;
    Vmin2 = max((Vini2*VolMax2-(FlowMax2+MaxSpill2)*TMax),VolMin2)/VolMax2;
    Vmax3 = min((Vini3*VolMax3+IT3(1)*TMax),VolMax3)/VolMax3;
    Vmin3 = max((Vini3*VolMax3-(FlowMax3+MaxSpill3)*TMax),VolMin3)/VolMax3;
    Vmax4 = min((Vini4*VolMax4+IT4(1)*TMax),VolMax4)/VolMax4;
    Vmin4 = max((Vini4*VolMax4-(FlowMax4+MaxSpill4)*TMax),VolMin4)/VolMax4;

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

    [vplot1,vplot3] = meshgrid(V1,V3); % For plotting surfaces.

    D = DMax*(0.5-(0.4)*sin(2*pi*time)); % In kW. (THIS (1))

    % D = DMax*(0.5)*time./time; % In kW. (THIS (2))

    % D = D01012017'*1000*0.8; % In kW. (THIS (3))
    % D = interp1(absTime,D,time); % Interpolation. (THIS (3))

    % ==================== Plot Demand (plot 6) ====================>>

    if ComputePlots == 1
        set(0,'CurrentFigure',1);
        hold on
        plot(time,D/DMax);
        xlabel('Normalized Time');
        ylabel('Normalized Demand');
        pause(0.001);
    elseif ComputePlots ~= 0
        disp('Choose ComputePlots between 0 or 1.');
        return;
    end

    % ==================== Matrices ====================>>

    u = zeros(length(V1),length(V3),length(V4));
    uV1 = zeros(length(V1),length(V3),length(V4));
    uV3 = zeros(length(V1),length(V3),length(V4));
    uV4 = zeros(length(V1),length(V3),length(V4));
    ControlsFMC = zeros(length(V1),length(V3),length(V4),9); % F,T1,S1,T2,S2,T3,T4,V2,V3.
    ControlsFMC_Aux = ControlsFMC;

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

    Khe1 = 0.3/10^6; % In M USD/MWh. Cost per energy of Bonete.
    Khe2 = 0/10^6; % In M USD/MWh.
    Khe3 = 0.3/10^6; % In M USD/MWh.
    Khe4 = 0.3/10^6; % In M USD/MWh.
    Kf = 130/10^6; % In M USD/MWh. Cost per energy of the fuel.

    CH1 = Khe1*eta*(H1(Vini1)-d1(1,0)-H2)/(3.6*10^6); % Water's value of Bonete.
    CH2 = Khe2*eta*(H2-d2(1,0)-H3(Vini3))/(3.6*10^6); % Water's value of Baygorria.
    CH3 = Khe3*eta*(H3(Vini3)-d3(1,0)-h03)/(3.6*10^6); % Water's value of Palmar.
    CH4 = Khe4*eta*(H4(Vini4)-d4(1,0)-h04)/(3.6*10^6); % Water's value of Salto Grande.
    CF = Kf/(3.6*10^6);

    % We divide over 3600 to pass from hours to seconds and over 1000 to move
    % from Mega to Kilo.

    olCF = CF*FuelMax;
    olCT1 = @(uV1,t) FlowMax1*(CH1 - Lambda21(t+Z21) - uV1/VolMax1 + Lambda21(t));
    olCS1 = @(uV1,t) MaxSpill1*(CH1 - Lambda21(t+Z21) - uV1/VolMax1 + Lambda21(t));
    olCT2 = @(uV2,t) FlowMax2*(CH2 - Lambda32(t+Z32) - uV2/VolMax2);
    olCS2 = @(uV2,t) MaxSpill2*(CH2 - Lambda32(t+Z32) - uV2/VolMax2);
    olCT3 = @(uV3,t) FlowMax3*(CH3 - uV3/VolMax3);
    olCT4 = @(uV4,t) FlowMax4*(CH4 - uV4/VolMax4);
    H0 = @(uV1,uV2,uV3,uV4,t) IT1(t)*uV1/VolMax1 + IT3(t)*uV3/VolMax3 + IT4(t)*uV4/VolMax4 - Lambda21(t)*IT2(t);
    olCV3 = @(uV3,t) MaxRel3*(Lambda32(t) + uV3/VolMax3);

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

                        if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                            uV1(v1,v3,v4) = (u(v1+1,v3,v4)-u(v1,v3,v4))/(dV1);
                        else
                            uV1(v1,v3,v4) = (u(v1,v3,v4)-u(v1-1,v3,v4))/(dV1);
                        end

                        if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,8) >= 0
                            uV3(v1,v3,v4) = (u(v1,v3+1,v4)-u(v1,v3,v4))/(dV3);
                        else
                            uV3(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3-1,v4))/(dV3);
                        end

                        if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0
                            uV4(v1,v3,v4) = (u(v1,v3,v4+1)-u(v1,v3,v4))/(dV4);
                        else
                            uV4(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3,v4-1))/(dV4);
                        end

                    elseif t < length(time) - 1

                        if v1 == 1
                            if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                                uV1(v1,v3,v4) = (u(v1+2,v3,v4)-u(v1+1,v3,v4))/(dV1);
                            else
                                uV1(v1,v3,v4) = (u(v1+1,v3,v4)-u(v1,v3,v4))/(dV1);
                            end
                        elseif v1 == length(V1)
                            if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                                uV1(v1,v3,v4) = (u(v1,v3,v4)-u(v1-1,v3,v4))/(dV1);
                            else
                                uV1(v1,v3,v4) = (u(v1-1,v3,v4)-u(v1-2,v3,v4))/(dV1);
                            end
                        elseif (v3==1 || v3==length(V3) || v4==1 || v4==length(V4)) && (v1~=1) && (v1~=length(V1))
                            if IT1(t+1)-FlowMax1*ControlsFMC_Aux(v1,v3,v4,2)-MaxSpill1*ControlsFMC_Aux(v1,v3,v4,3) >= 0
                                uV1(v1,v3,v4) = (u(v1+1,v3,v4)-u(v1,v3,v4))/(dV1);
                            else
                                uV1(v1,v3,v4) = (u(v1,v3,v4)-u(v1-1,v3,v4))/(dV1);
                            end
                        end

                        if v3 == 1
                            if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,8) >= 0
                                uV3(v1,v3,v4) = (u(v1,v3+2,v4)-u(v1,v3+1,v4))/(dV3);
                            else
                                uV3(v1,v3,v4) = (u(v1,v3+1,v4)-u(v1,v3,v4))/(dV3);
                            end
                        elseif v3 == length(V3)
                            if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,8) >= 0
                                uV3(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3-1,v4))/(dV3);
                            else
                                uV3(v1,v3,v4) = (u(v1,v3-1,v4)-u(v1,v3-2,v4))/(dV3);
                            end
                        elseif (v1==1 || v1==length(V1) || v4==1 || v4==length(V4)) && (v3~=1) && (v3~=length(V3))
                            if IT3(t+1)-FlowMax3*ControlsFMC_Aux(v1,v3,v4,6)+MaxV3*ControlsFMC_Aux(v1,v3,v4,8) >= 0
                                uV3(v1,v3,v4) = (u(v1,v3+1,v4)-u(v1,v3,v4))/(dV3);
                            else
                                uV3(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3-1,v4))/(dV3);
                            end
                        end

                        if v4 == 1
                            if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0
                                uV4(v1,v3,v4) = (u(v1,v3,v4+2)-u(v1,v3,v4+1))/(dV4);
                            else
                                uV4(v1,v3,v4) = (u(v1,v3,v4+1)-u(v1,v3,v4))/(dV4);
                            end
                        elseif v4 == length(V4)
                            if IT4(t+1)-FlowMax4*ControlsFMC_Aux(v1,v3,v4,7) >= 0
                                uV4(v1,v3,v4) = (u(v1,v3,v4)-u(v1,v3,v4-1))/(dV4);
                            else
                                uV4(v1,v3,v4) = (u(v1,v3,v4-1)-u(v1,v3,v4-2))/(dV4);
                            end
                        elseif (v1==1 || v1==length(V1) || v3==1 || v3==length(V3)) && (v4~=1) && (v4~=length(V4))
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

        % ==================== Controls Section ====================>>

        ControlsFMC_Aux = ControlsFMC;

        for v1 = 1:length(V1)

            for v3 = 1:length(V3)

                for v4 = 1:length(V4)

                    ControlsFMC(v1,v3,v4,1:8) = Quad_FMC_Erik(Q,b(V1(v1),V3(v3),V4(v4)),c(t),...
                        d(uV1(v1,v3,v4),0,uV3(v1,v3,v4),uV4(v1,v3,v4),t),k(t),options); % dV2 = 0.

                    if olCV3(uV3(v1,v3,v4),t) >= 0
                        ControlsFMC(v1,v3,v4,9) = 0;
                    else
                        ControlsFMC(v1,v3,v4,9) = 0; % This is 1.
                    end

                end

            end

        end

        % ==================== Time Step ====================>>

        if t ~= length(time)

            for v1 = 1:length(V1)

                for v3 = 1:length(V3)

                    for v4 = 1:length(V4)

                        Controls(1:9) = ControlsFMC_Aux(v1,v3,v4,:);

                        u(v1,v3,v4) = u(v1,v3,v4) + dt*TMax*...
                            (Controls*...
                            [d(uV1(v1,v3,v4),0,uV3(v1,v3,v4),uV4(v1,v3,v4),t);0;olCV3(uV3(v1,v3,v4),t)]+...
                            H0(uV1(v1,v3,v4),0,uV3(v1,v3,v4),uV4(v1,v3,v4),t));

                    end

                end

            end

        end

        % ==================== Plotting over time ====================>>

        if ComputePlots == 1

            set(0,'CurrentFigure',1);
            vline(time(t),'r',' ');
            title(['Demand (time=',num2str(t),')']);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Demand_',num2str(expansion)],'epsc');
            end

            set(0,'CurrentFigure',2);
            clf(2);
            hold on;
            for index=1:floor(length(V4)/2):length(V4)
                    S = surf(vplot1,vplot3,squeeze(u(:,:,index))');
                    S.EdgeColor = 'none';
            end
            xlabel('Vol. Bonete');
            ylabel('Vol. Palmar');
            zlabel('M USD');
            title(['Optimal cost function (time=',num2str(t),')']);
            xlim([Vmin1 Vmax1]);
            ylim([Vmin3 Vmax3]);
            zlim auto; box; view(140,20);
            if SavePlots == 1 && t == 1
                saveas(gcf,['Cost_',num2str(expansion)],'epsc');
            end

        elseif ComputePlots ~= 0
            disp('Choose ComputePlots between 0 or 1.');
            return;
        end    
        
        fprintf('Discretizations = %d. Completed : %d',NV4,100*(1-t/length(time)));
        fprintf('\n');
        pause(0.1);

    end

    CP = floor(length(V1)/2) + 1;
    Uini = u(CP,CP,CP);

end