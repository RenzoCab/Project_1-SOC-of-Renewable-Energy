function [] = Real_Data(SaveFig,Plot,SaveDay)

    Data = xlsread([pwd '/Data/Real_January-February_2019.xlsx'],1);
    Fuel = xlsread([pwd '/Data/Fuel_From_SimSEE.xlsx'],1);
    Agua_Vals = xlsread([pwd '/Data/ValsAguaProm_SimSEE_Enero_Febrero_2019.xlsx'],1);
    SG = xlsread([pwd '/Data/Salto_Grande.xlsx'],1);
    Bonete = xlsread([pwd '/Data/Represas_Bonete.xlsx'],1);
    Baygorria = xlsread([pwd '/Data/Represas_Baygorria.xlsx'],1);
    Palmar = xlsread([pwd '/Data/Represas_Palmar.xlsx'],1);
    All = xlsread([pwd '/Data/All_Dams_JenFeb.xlsx'],1);

    F = {'Fecha' 'Salto Grande' 'Bonete' 'Baygorria' 'Palmar' 'Eólica' 'Solar' 'Térmica' 'Biomasa' 'Imp.Arg' 'Imp.Br.Riv' ...
        'Imp.Br.Mel' 'Demanda' 'Exp_Intercon_BR_MELO' 'Exp_Intercon_BR_RIVERA' 'Exp_Intercon_ARG' 'Imp_Intercon_BR_MELO' 'Imp_Intercon_BR_RIVERA' 'Imp_Intercon_AG_Imp'};

    Fecha = {'01/01/2019','02/01/2019','03/01/2019','04/01/2019','05/01/2019','06/01/2019','07/01/2019','08/01/2019','09/01/2019','10/01/2019',...
        '11/01/2019','12/01/2019','13/01/2019','14/01/2019','15/01/2019','16/01/2019','17/01/2019','18/01/2019','19/01/2019','20/01/2019',...
        '21/01/2019','22/01/2019','23/01/2019','24/01/2019','25/01/2019','26/01/2019','27/01/2019','28/01/2019','29/01/2019','30/01/2019',...
        '31/01/2019','01/02/2019','02/02/2019','03/02/2019','04/02/2019','05/02/2019','06/02/2019','07/02/2019','08/02/2019','09/02/2019','10/02/2019',...
        '11/02/2019','12/02/2019','13/02/2019','14/02/2019','15/02/2019','16/02/2019','17/02/2019','18/02/2019','19/02/2019','20/02/2019',...
        '21/02/2019','22/02/2019','23/02/2019','24/02/2019','25/02/2019','26/02/2019','27/02/2019','28/02/2019'};

    T = (length(Data(:,1))-1-60)/6/24;
    dt = 1/24/6;
    time = 0:dt:1;

    for i = 0:T-1
        if Plot == 1
            figure('Position', [10 10 900 420]);
            hh = area(time,[Data(60+1+145*i-i:60+145+145*i-i,9),Data(60+1+145*i-i:60+145+145*i-i,7),Data(60+1+145*i-i:60+145+145*i-i,6),...
                Data(60+1+145*i-i:60+145+145*i-i,2),Data(60+1+145*i-i:60+145+145*i-i,3),Data(60+1+145*i-i:60+145+145*i-i,4),Data(60+1+145*i-i:60+145+145*i-i,5),Data(60+1+145*i-i:60+145+145*i-i,8)]);
            hold on;
            P = plot(time,Data(60+1+145*i-i:60+145+145*i-i,13),'k'); P.LineWidth = 1;
            P = plot(time,Data(60+1+145*i-i:60+145+145*i-i,13)+Data(60+1+145*i-i:60+145+145*i-i,14)+Data(60+1+145*i-i:60+145+145*i-i,15)+Data(60+1+145*i-i:60+145+145*i-i,16),'m');
            P.LineWidth = 1;
            grid minor;
            legend(F{9},F{7},F{6},F{2:5},F{8},F{13},'Demanda + Exp.','location','eastoutside');
            hh(1).FaceColor = [1 .5 0]; % Biomasa.
            hh(2).FaceColor = [1 1 .5]; % Solar.
            hh(3).FaceColor = [0 1 0]; % Wind.
            hh(4).FaceColor = [.5 .5 1]; % Bonete.
            hh(5).FaceColor = [.5 .8 1]; % Baygorria.
            hh(6).FaceColor = [.2 .2 1]; % Palmar.
            hh(7).FaceColor = [.8 .8 1]; % SG.
            hh(8).FaceColor = [1 0 0]; % Fuel.
            ylabel('MW');
            xlabel('Time');
            title(['Historical Production ',Fecha(i+1)]);
        end
        if 0 == exist('Historical','dir')
            mkdir('Historical');
        end
        if SaveDay == 1
            Inflow2 = Baygorria(9,T-i) - (1/4)*(Bonete(5,T-i+1)+Bonete(6,T-i+1)) - (3/4)*(Bonete(5,T-i)+Bonete(6,T-i));
            Inflow3 = Palmar(9,T-i) - (10/24)*(Baygorria(5,T-i+1)+Baygorria(6,T-i+1)) - (14/24)*(Baygorria(5,T-i)+Baygorria(6,T-i));
            Matrix = {Data(60+1+145*i-i:60+145+145*i-i,:),Fuel(i+1,2:9),Agua_Vals(i+1,2:4),...
                [Bonete(2,T-i),Bonete(9,T-i)],[Baygorria(2,T-i),Inflow2],[Palmar(2,T-i),Inflow3],...
                [SG(i+2,3),SG(i+1,4)],[All(1+24*i:25+24*i,8),All(1+24*i:25+24*i,11)...
                ,All(1+24*i:25+24*i,13:16),ones(25,1)*mean(All(1+24*i:25+24*i,17)),ones(25,1)*mean(All(1+24*i:25+24*i,18))],...
                Data(1+145*i-i:1+60+145*i-i,3:4)};
            save([pwd '/Historical/Day_',num2str(i+1),'_2019.mat'],'Matrix');
        end
        if SaveFig == 1
            saveas(gcf,[pwd '/Historical/',num2str(i+1)],'epsc');
        end
    end

end