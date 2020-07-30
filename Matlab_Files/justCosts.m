close all;
clear all;
clc;

Day = '20180302';
finalDate = '20190101';
set(0,'DefaultFigureVisible','off');
figure(1);

if 0 == exist([pwd '/Simulations/plotsCosts'],'dir')
    mkdir([pwd '/Simulations/plotsCosts']);
end

while not(strcmp(Day,finalDate))
    
    time = [0,1];
    Legend = {'Bonete','Baygorria','Palmar','Salto Grande','Motores Batlle','PTA','PTB','CTR'};
    Colors = {[.5 .8 1],[.2 .2 1],[.8 .8 1],[.5 .5 1],[1 .6 .6],[1 .3 .5],[145/255 0 211/255],[1 0 0]};
    load([pwd '/Historical/dailyData/Day_',Day,'.mat']);
    Costs(5:8) = [Matrix{2}(3),Matrix{2}(5),Matrix{2}(7),Matrix{2}(1)]; % USD/MWh.
    Costs(1) = Matrix{3}(3); % In USD/MWh. Cost per energy of Bonete.
    Costs(2) = Matrix{3}(4); % In USD/MWh. Cost per energy of Baygorria.
    Costs(3) = Matrix{3}(2); % In USD/MWh. Cost per energy of Palmar.
    Costs(4) = Matrix{3}(1); % In USD/MWh. Cost per energy of SG.
    
%     hh(1).FaceColor = [1 .5 0]; % Biomasa.
%     hh(2).FaceColor = [1 1 .5]; % Solar.
%     hh(3).FaceColor = [0 1 0]; % Wind.
%     hh(4).FaceColor = [.5 .5 1]; % SG.
%     hh(5).FaceColor = [.5 .8 1]; % Bonete.
%     hh(6).FaceColor = [.2 .2 1]; % Baygorria.
%     hh(7).FaceColor = [.8 .8 1]; % Palmar.
%     hh(8).FaceColor = [1 .6 .6]; % T1.
%     hh(9).FaceColor = [1 .3 .5]; % T2.
%     hh(10).FaceColor = [145/255 0 211/255;]; % T3.
%     hh(11).FaceColor = [1 0 0]; % T4.
%     hh(12).FaceColor = [0 0 0]; % T5.
    
    set(0,'CurrentFigure',1); clf(1);
    hold on;
    for cont = 1:length(Costs)
        P = plot(time,ones(1,length(time))*Costs(cont));
        P.LineWidth = 3;
        P.Color = Colors{cont};
    end
    grid minor;
    title('Daily Water and Fossil Fuels Costs');
    ylabel('USD/MWh');
    legend(Legend);
    set(gca,'Box','on');
    
	saveas(gcf,[pwd '/Simulations/plotsCosts/costs_',Day],'epsc');

    Day = nextDay(Day);
    
end

set(0,'DefaultFigureVisible','on');