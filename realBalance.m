%% Function realBalance:
% This function creates the plots with the real historical balance. It uses
% as an input a date in format yyyymmdd (string) and does not has output.
% However, it displays the plot with the real power balance and saves it.

function [] = realBalancePlot(dateString)

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
    
    data_1_1 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos/',dateString,'.ods'],'GPF',''));
    data_1_2 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos/',dateString,'.ods'],'Eólica',''));
    data_1_3 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos/',dateString,'.ods'],'Solar',''));
    data_1_4 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos/',dateString,'.ods'],'Térmica',''));
    data_1_5 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos/',dateString,'.ods'],'Biomasa',''));
    data_1_6 = cell2mat(loadods(['../Python/Represas_Data_2/ADME_Historicos/',dateString,'.ods'],'Intercambios.',''));
    
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
    ylabel('MW');
    xlabel('Time');
    title(['Daily Historical Production ',dateString_2,' (ADME Data)']);

end