%% Function realBalancePlotCont:
% This function creates the plots with the real historical balance, but 
% only showing the controllable sources. It uses
% as an input a date in format yyyymmdd (string) and does not has output.
% However, it displays the plot with the real power balance and saves it.

function [] = realBalancePlotCont(dateString)

    try
        dateString;
        datetime(dateString,'InputFormat','yyyyMMdd');
    catch
        warning('The input dateString is not specified, or it has an incorrect format. The format is yyyymmdd string, i.e., 20190101.');
        warning('By default the date 20190101 will be use (this is useful for debugging).');
        dateString = '20190101';
    end
    
    dateTime = datetime(dateString,'InputFormat','yyyyMMdd');
    dateString_2 = datestr(dateTime,'dd/mm/yyyy');
    
    dt = 1/24/6;
    time = 0:dt:1;
    Termicas = {'Motores Batlle (Fossil)','PTA8 (Fossil)','PTB (Fossil)','CTR (Fossil)'};
    Dams = {'Salto Grande (Hydro)','Bonete (Hydro)','Baygorria (Hydro)','Palmar (Hydro)'};
    
    data_1_1 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString,'.ods'],'GPF',''));
    data_1_2 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString,'.ods'],'Eolica',''));
    data_1_3 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString,'.ods'],'Solar',''));
    data_1_4 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString,'.ods'],'Termica',''));
    data_1_5 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString,'.ods'],'Biomasa',''));
    data_1_6 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString,'.ods'],'Intercambios.',''));
    
    data_1_1(end,:) = data_1_1(end-1,:);
    data_1_4(end,:) = data_1_4(end-1,:);
    data_1_6(end,:) = data_1_6(end-1,:);
    % I do this only to avoid loading the next-day file, which has
    % the information correctly at time 00:00.
    
    Zenda = data_1_4(:,1);
    PTA6 = data_1_4(:,2);
    Motores = data_1_4(:,3);
    CTR = data_1_4(:,4);
    PTB = data_1_4(:,5);
    PTA78 = data_1_4(:,5);
    PTA8 = PTA6 + PTA78;
    
	figure('Position', [10 10 900 420]);
    
    plotArea = area(time,[data_1_1(:,1),data_1_1(:,2),data_1_1(:,3),data_1_1(:,4),Motores,PTA8,PTB,CTR]);
    hold on;
    P = plot(time,data_1_1(:,12)+data_1_6(:,1)+data_1_6(:,2)+data_1_6(:,3)...
        -data_1_1(:,5)-data_1_1(:,6)-data_1_1(:,8),'m');
    P.LineWidth = 1;
    grid minor;
    legend(Dams{:},Termicas{:},'Effective Demand','location','eastoutside');
    plotArea(1).FaceColor = [.5 .5 1]; % SG.
    plotArea(2).FaceColor = [.5 .8 1]; % Bonete.
    plotArea(3).FaceColor = [.2 .2 1]; % Baygorria.
    plotArea(4).FaceColor = [.8 .8 1]; % Palmar.
    plotArea(5).FaceColor = [1 .6 .6]; % T1 (Motores).
    plotArea(6).FaceColor = [1 .3 .5]; % T2 (PTA8).
    plotArea(7).FaceColor = [145/255 0 211/255;]; % T3 (PTB).
    plotArea(8).FaceColor = [1 0 0]; % T4 (CTR).
    ylabel('Power (MW)');
    xlabel('Time');
    title(['Daily Historical Controllable Production ',dateString_2,' (ADME Data)']);
    
    if 0 == exist('Historical','dir')
    	mkdir('Historical');
    end
    if 0 == exist('Historical/realBalancePlot','dir')
    	mkdir('Historical/realBalancePlot');
    end
    
    saveas(gcf,[pwd '/Historical/realBalancePlot/contOnly_',dateString],'epsc');

end