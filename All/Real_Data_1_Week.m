function [] = Real_Data_1_Week(SaveFig,Plot,SaveDay)

    Data = xlsread([pwd '/Data/Real_January-February_2019.xlsx'],1);
    Fuel = xlsread([pwd '/Data/Fuel_From_SimSEE.xlsx'],1);
    Agua_Vals = xlsread([pwd '/Data/ValsAguaProm_SimSEE_Enero_Febrero_2019.xlsx'],1);
    SG = xlsread([pwd '/Data/Salto_Grande.xlsx'],1);
    Bonete = xlsread([pwd '/Data/Represas_Bonete.xlsx'],1);
    Baygorria = xlsread([pwd '/Data/Represas_Baygorria.xlsx'],1);
    Palmar = xlsread([pwd '/Data/Represas_Palmar.xlsx'],1);
    All = xlsread([pwd '/Data/All_Dams_JenFeb.xlsx'],1);

    F = {'Fecha' 'Salto Grande' 'Bonete' 'Baygorria' 'Palmar' 'Wind Power' 'Solar Power' 'Fuel' 'Biomass' 'Imp.Arg' 'Imp.Br.Riv' ...
        'Imp.Br.Mel' 'Demand' 'Exp_Intercon_BR_MELO' 'Exp_Intercon_BR_RIVERA' 'Exp_Intercon_ARG' 'Imp_Intercon_BR_MELO' 'Imp_Intercon_BR_RIVERA' 'Imp_Intercon_AG_Imp'};

    Fecha = {'01/01/2019','02/01/2019','03/01/2019','04/01/2019','05/01/2019','06/01/2019','07/01/2019','08/01/2019','09/01/2019','10/01/2019',...
        '11/01/2019','12/01/2019','13/01/2019','14/01/2019','15/01/2019','16/01/2019','17/01/2019','18/01/2019','19/01/2019','20/01/2019',...
        '21/01/2019','22/01/2019','23/01/2019','24/01/2019','25/01/2019','26/01/2019','27/01/2019','28/01/2019','29/01/2019','30/01/2019',...
        '31/01/2019','01/02/2019','02/02/2019','03/02/2019','04/02/2019','05/02/2019','06/02/2019','07/02/2019','08/02/2019','09/02/2019','10/02/2019',...
        '11/02/2019','12/02/2019','13/02/2019','14/02/2019','15/02/2019','16/02/2019','17/02/2019','18/02/2019','19/02/2019','20/02/2019',...
        '21/02/2019','22/02/2019','23/02/2019','24/02/2019','25/02/2019','26/02/2019','27/02/2019','28/02/2019'};

    T = (length(Data(:,1))-1-60)/6/24;
    dt = 1/24/6;
    time = 0:dt:7;

    Termicas = {'Motores Batlle','PTA','PTB','CTR'};
    Dams = {'Salto Grande','Bonete','Baygorria','Palmar'};
    
    if Plot == 1
        for i = 0:0
            figure('Position', [10 10 1500 420]);
            CTR = interp1(linspace(1,24*7+1,24*7+1),All(1:24*7+1,13),linspace(1,24*7+1,144*7+1));
            PTA = interp1(linspace(1,24*7+1,24*7+1),All(1:24*7+1,14),linspace(1,24*7+1,144*7+1));
            PTB = interp1(linspace(1,24*7+1,24*7+1),All(1:24*7+1,15),linspace(1,24*7+1,144*7+1));
            Motores_Batlle = interp1(linspace(1,24*7+1,24*7+1),All(1:24*7+1,13),linspace(1,24*7+1,144*7+1));
            hh = area(time,[Data(60+1:60+144*7+1,9),Data(60+1:60+144*7+1,7),Data(60+1:60+144*7+1,6),...
            Data(60+1:60+144*7+1,2),Data(60+1:60+144*7+1,3),Data(60+1:60+144*7+1,4),Data(60+1:60+144*7+1,5),Motores_Batlle',...
            PTA',PTB',CTR']);
            hold on;
            P = plot(time,Data(60+1:60+144*7+1,13),'k'); P.LineWidth = 1;
            P = plot(time,Data(60+1:60+144*7+1,13)+Data(60+1:60+144*7+1,14)+Data(60+1:60+144*7+1,15)+Data(60+1:60+144*7+1,16),'m');
            P.LineWidth = 1;
            grid minor;
            legend(F{9},F{7},F{6},F{2:5},Termicas{:},F{13},'Demand + Exp','location','eastoutside');
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
            ylabel('Power (MW)');
            xlabel('Time');
            title(['Historical Production ',Fecha(i+1)]);
            if 0 == exist('Historical_1W','dir')
                mkdir('Historical_1W');
            end
            if SaveDay == 1
                Inflow2 = Baygorria(9,T-6:T) - (1/4)*(Bonete(5,T-5:T+1)+Bonete(6,T-5:T+1)) - (3/4)*(Bonete(5,T-6:T)+Bonete(6,T-6:T));
                Inflow3 = Palmar(9,T-6:T) - (10/24)*(Baygorria(5,T-5:T+1)+Baygorria(6,T-5:T+1)) - (14/24)*(Baygorria(5,T-6:T)+Baygorria(6,T-6:T));
                Matrix = {Data(60+1:60+144*7+1,:),Fuel(1:7,2:9),Agua_Vals(1:7,2:4),...
                    [Bonete(2,T-6:T)',Bonete(9,T-6:T)'],[Baygorria(2,T-6:T)',Inflow2'],[Palmar(2,T-6:T)',Inflow3'],...
                    [SG(2:8,3),SG(1:7,4)],[All(1:24*7+1,8),All(1:24*7+1,11)...
                    ,All(1:24*7+1,13:16),ones(24*7+1,1)*mean(All(1:24*7+1,17)),ones(24*7+1,1)*mean(All(1:24*7+1,18))],...
                    Data(1:1+60,3:4)};
                save([pwd '/Historical_1W/Day_',num2str(i+1),'_2019.mat'],'Matrix');
            end
            if SaveFig == 1
                saveas(gcf,[pwd '/Historical_1W/A_',num2str(i+1)],'epsc');
            end
        end
        close all;
        pause(0.1);
        for i = 0:0
            figure('Position', [10 10 900 420]);
            CTR = interp1(linspace(1,25,25),All(1+24*i:25+24*i,13),linspace(1,25,145));
            PTA = interp1(linspace(1,25,25),All(1+24*i:25+24*i,14),linspace(1,25,145));
            PTB = interp1(linspace(1,25,25),All(1+24*i:25+24*i,15),linspace(1,25,145));
            Motores_Batlle = interp1(linspace(1,25,25),All(1+24*i:25+24*i,13),linspace(1,25,145));
            hh = area(time,[Data(60+1+145*i-i:60+145+145*i-i,2),Data(60+1+145*i-i:60+145+145*i-i,3),...
                Data(60+1+145*i-i:60+145+145*i-i,4),Data(60+1+145*i-i:60+145+145*i-i,5),Motores_Batlle',PTA',PTB',CTR']);
            hold on;
            P = plot(time,Data(60+1+145*i-i:60+145+145*i-i,13)+Data(60+1+145*i-i:60+145+145*i-i,14)+Data(60+1+145*i-i:60+145+145*i-i,15)...
                +Data(60+1+145*i-i:60+145+145*i-i,16)-Data(60+1+145*i-i:60+145+145*i-i,9)-Data(60+1+145*i-i:60+145+145*i-i,7)-Data(60+1+145*i-i:60+145+145*i-i,6),'m');
            P.LineWidth = 1;
            grid minor;
            legend(Dams{:},Termicas{:},'Demand + Exp','location','eastoutside');
            hh(1).FaceColor = [.5 .5 1]; % Bonete.
            hh(2).FaceColor = [.5 .8 1]; % Baygorria.
            hh(3).FaceColor = [.2 .2 1]; % Palmar.
            hh(4).FaceColor = [.8 .8 1]; % SG.
            hh(5).FaceColor = [1 .6 .6]; % T1.
            hh(6).FaceColor = [1 .3 .5]; % T2.
            hh(7).FaceColor = [145/255 0 211/255;]; % T3.
            hh(8).FaceColor = [1 0 0]; % T4.
            ylabel('MW');
            xlabel('Time');
            title(['Historical Production ',Fecha(i+1)]);
            if SaveFig == 1
                saveas(gcf,[pwd '/Historical_1W/B_',num2str(i+1)],'epsc');
            end
        end
        close all;
        pause(0.1);
        for i = 0:0
            figure('Position', [10 10 900 420]);
            hh = area(time,[Data(60+1+145*i-i:60+145+145*i-i,9),Data(60+1+145*i-i:60+145+145*i-i,7),Data(60+1+145*i-i:60+145+145*i-i,6),...
                Data(60+1+145*i-i:60+145+145*i-i,2),Data(60+1+145*i-i:60+145+145*i-i,3),Data(60+1+145*i-i:60+145+145*i-i,4),Data(60+1+145*i-i:60+145+145*i-i,5),Data(60+1+145*i-i:60+145+145*i-i,8)]);
            hold on;
            P = plot(time,Data(60+1+145*i-i:60+145+145*i-i,13),'k'); P.LineWidth = 1;
            P = plot(time,Data(60+1+145*i-i:60+145+145*i-i,13)+Data(60+1+145*i-i:60+145+145*i-i,14)+Data(60+1+145*i-i:60+145+145*i-i,15)+Data(60+1+145*i-i:60+145+145*i-i,16),'m');
            P.LineWidth = 1;
            grid minor;
            legend(F{9},F{7},F{6},F{2:5},F{8},F{13},'Demand + Exp','location','eastoutside');
            hh(1).FaceColor = [1 .5 0]; % Biomasa.
            hh(2).FaceColor = [1 1 .5]; % Solar.
            hh(3).FaceColor = [0 1 0]; % Wind.
            hh(4).FaceColor = [.5 .5 1]; % Bonete.
            hh(5).FaceColor = [.5 .8 1]; % Baygorria.
            hh(6).FaceColor = [.2 .2 1]; % Palmar.
            hh(7).FaceColor = [.8 .8 1]; % SG.
            hh(8).FaceColor = [1 0 0]; % Fuel.
            ylabel('Power (MW)');
            xlabel('Time');
            title(['Historical Production ',Fecha(i+1)]);
            if SaveFig == 1
                saveas(gcf,[pwd '/Historical_1W/C_',num2str(i+1)],'epsc');
            end
        end
        close all;
        pause(0.1);
        for i = 0:0
            figure('Position', [10 10 900 420]);
            hh = area(time,[Data(60+1+145*i-i:60+145+145*i-i,2),Data(60+1+145*i-i:60+145+145*i-i,3),...
                Data(60+1+145*i-i:60+145+145*i-i,4),Data(60+1+145*i-i:60+145+145*i-i,5),Data(60+1+145*i-i:60+145+145*i-i,8)]);
            hold on;
            P = plot(time,Data(60+1+145*i-i:60+145+145*i-i,13)+Data(60+1+145*i-i:60+145+145*i-i,14)+Data(60+1+145*i-i:60+145+145*i-i,15)...
                +Data(60+1+145*i-i:60+145+145*i-i,16)-Data(60+1+145*i-i:60+145+145*i-i,9)-Data(60+1+145*i-i:60+145+145*i-i,7)-Data(60+1+145*i-i:60+145+145*i-i,6),'m');
            P.LineWidth = 1;
            grid minor;
            legend(Dams{:},F{8},'Demand + Exp','location','eastoutside');
            hh(1).FaceColor = [.5 .5 1]; % Bonete.
            hh(2).FaceColor = [.5 .8 1]; % Baygorria.
            hh(3).FaceColor = [.2 .2 1]; % Palmar.
            hh(4).FaceColor = [.8 .8 1]; % SG.
            hh(5).FaceColor = [1 0 0]; % Fuel.
            ylabel('MW');
            xlabel('Time');
            title(['Historical Production ',Fecha(i+1)]);
            if SaveFig == 1
                saveas(gcf,[pwd '/Historical_1W/D_',num2str(i+1)],'epsc');
            end
        end
    end
    
    %% For computing the total energy:
    
    for i = 0:T-1
        
        Exp = 8;
        
        % Energy: Bonete, Baygorria, Palmar, SG, T1..4.
        CTR = interp1(linspace(1,25,25),All(1+24*i:25+24*i,13),linspace(1,25,145));
        PTA = interp1(linspace(1,25,25),All(1+24*i:25+24*i,14),linspace(1,25,145));
        PTB = interp1(linspace(1,25,25),All(1+24*i:25+24*i,15),linspace(1,25,145));
        Motores_Batlle = interp1(linspace(1,25,25),All(1+24*i:25+24*i,13),linspace(1,25,145));
        Total_Energy = {sum(Data(60+1+145*i-i:60+145+145*i-i,3))/6,sum(Data(60+1+145*i-i:60+145+145*i-i,4))/6,...
            sum(Data(60+1+145*i-i:60+145+145*i-i,5))/6,sum(Data(60+1+145*i-i:60+145+145*i-i,2))/6,...
            sum(Motores_Batlle)/6,sum(PTA)/6,sum(PTB)/6,sum(CTR)/6};
        
        % Controls:
        Control_Bonete = interp1(linspace(1,25,145),Data(60+1+145*i-i:60+145+145*i-i,3),linspace(1,25,2^Exp+1));
        Control_Baygorria = interp1(linspace(1,25,145),Data(60+1+145*i-i:60+145+145*i-i,3),linspace(1,25,2^Exp+1));
        Control_Palmar = interp1(linspace(1,25,145),Data(60+1+145*i-i:60+145+145*i-i,5),linspace(1,25,2^Exp+1));
        Control_SG = interp1(linspace(1,25,145),Data(60+1+145*i-i:60+145+145*i-i,2),linspace(1,25,2^Exp+1));
        Control_Motores_Batlle = interp1(linspace(1,25,145),Motores_Batlle,linspace(1,25,2^Exp+1));
        Control_PTA = interp1(linspace(1,25,145),PTA,linspace(1,25,2^Exp+1));
        Control_PTB = interp1(linspace(1,25,145),PTB,linspace(1,25,2^Exp+1));
        Control_CTR = interp1(linspace(1,25,145),CTR,linspace(1,25,2^Exp+1));
        Controls = {Control_Bonete,Control_Baygorria,Control_Palmar,Control_SG,...
            Control_Motores_Batlle,Control_PTA,Control_PTB,Control_CTR};
        
        % To save:
        if SaveDay == 1
            Save = {Total_Energy,Controls};
            save([pwd '/Historical_1W/Energy_Controls_',num2str(i+1),'_2019.mat'],'Save');
        end
        
    end
    
end