function [] = New_Real_Data(DateVector)

% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email: renzo.caballerorosas@kaust.edu.sa caballerorenzo@hotmail.com
% Website: None.
% June 2019; Last revision: 16/12/2019.

% DateVector must have the format DateVector = [yyyy,mm,dd], we do not care
% about the time as every simulation starts at 00:00. The output is a
% matrix (called Matrix) which is saved in a file yyyymmdd.mat.

    %%%%%%%%%% Data #1: Los nombres de las hojas tienen la informacion necesaria.

    % ==================
    %DateVector = [DateVector,0,0,0]; % Here we add the time 00:00:00 if DateVector has format [yyyy,mm,dd].
    %DateVector = [2019,4,2,0,0,0]; % 02/04/2019. If we want to do some test, we can uncomment this line.
    dateFormat = datestr(DateVector);
    dateString_1 = datestr(dateFormat,'yyyymmdd');
    % dateString is equal to '20190402'.
    dateString_2 = datestr(dateFormat,'yyyymmdd HH:MM:SS.FFF');
    % dataString_2 is equal to '20190402 00:00:00.000'.
    dateString_3 = datestr(dateFormat,'dd/mm/yyyy');
    % dataString_3 is equal to '02/04/2019'.
    dateString_4 = datestr(dateFormat,'yymmdd');
    % dateString is equal to '190402'.
    % ==================
        
    data_1_1 = loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString_1,'.ods'],'GPF','');
    data_1_2 = loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString_1,'.ods'],'Eolica','');
    data_1_3 = loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString_1,'.ods'],'Solar','');
    data_1_4 = loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString_1,'.ods'],'Termica','');
    data_1_5 = loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString_1,'.ods'],'Biomasa','');
    data_1_6 = loadods(['../Python/Represas_Data_2/ADME_Historicos_Corrected/',dateString_1,'.ods'],'Intercambios.','');
    
    data_1 = [cell2mat(data_1_1),cell2mat(data_1_2),cell2mat(data_1_3),cell2mat(data_1_4),...
        cell2mat(data_1_5),cell2mat(data_1_6)]; % To verify!!!
    data_1(end,:) = data_1(end-1,:);
    
    %%%%%%%%%% Data #2:
    
    % ==================
    excelVar = 1; % In Excel, '1' corresponds to 1 of January 1900.
    matlabVar = x2mdate(excelVar); % From Excel date-number to Matlab date-number.
    excelVar = m2xdate(matlabVar); % From Matlab date-number to Excel date-number.
    matlabDate = datestr(matlabVar); % From Matlab date-number to Matlab date-string.
    matlabDate = datenum(matlabDate); % From Matlab date-string to Matlab date-number.
    % ==================
    
    names_ADME = {'CC540','CTR','Motores','PTA_78','PTI_GO'};
    excelDate = m2xdate(datenum(dateFormat)); % Date in the numeric format of Excel.
    data_2_1_1 = xlsread(['./Data_Valores/Termicas_ADME.xlsx'],1); % Ciclo Combinado (PTB).
    data_2_1_2 = xlsread(['./Data_Valores/Termicas_ADME.xlsx'],2); % CTR.
    data_2_1_3 = xlsread(['./Data_Valores/Termicas_ADME.xlsx'],3); % Motores Central Batlle.
    data_2_1_4 = xlsread(['./Data_Valores/Termicas_ADME.xlsx'],4); % PTA 7 y 8.
    data_2_1_5 = xlsread(['./Data_Valores/Termicas_ADME.xlsx'],5); % Punta del Tigre (PTA6).
    
    CTR_Cost = data_2_1_2(data_2_1_2==excelDate,2);
    Mot_Cost = data_2_1_3(data_2_1_3==excelDate,2);
    PTA6_Cost = data_2_1_5(data_2_1_5==excelDate,2);
    PTB_Cost = data_2_1_1(data_2_1_1==excelDate,2);
    PTA78_Cost = data_2_1_4(data_2_1_4==excelDate,2);
    
    % This section is not being used.
    % PTA8 = PTA6 + PTA7&8.
    % Potencias maximas:
    %PTA6 = 300;
    %PTA78 = 50;
    %PTA8 = PTA6 + PTA78;
    %PTB = 550;
    %MB = 80;
    %CTR = 200;
    
    date_time = datetime(dateString_1,'InputFormat','yyyyMMdd'); % This is the date in datetime format.
    date_time_before = date_time - days(1);
    dateString_3_before = datestr(date_time_before,'dd/mm/yyyy');
    dateString_1_before = datestr(date_time_before,'yyyymmdd');
    
    data_2_2_1 = readtable(['./Data_Valores/Termicas_SimSEE.xlsx'],'Sheet',1); % Potencia maxima por turbina.
    aux_1_1 = data_2_2_1.Var1;
    aux_1_2 = data_2_2_1.Fecha;
    for i = 1:length(names_ADME)
        index_1{i} = find(strcmp(aux_1_1,names_ADME{i}));
        dates_1{i} = aux_1_2(index_1{i});
        aux = {};
        for j = 1:length(dates_1{i})
            aux{j} = datetime(dates_1{i}{j},'InputFormat','dd-MM-yyyy');
        end
        dates_time_1{i} = aux;
        
        for j = 1:length(dates_time_1{i})
            if date_time >= dates_time_1{i}{j}
                maxTurbine(i) = data_2_2_1.Power(index_1{i}(j));
            end
        end
    end
    
    data_2_2_2 = readtable(['./Data_Valores/Termicas_SimSEE.xlsx'],'Sheet',2); % Turbinas disponibles y mantenimiento programado.
    aux_2_1 = data_2_2_2.Var1;
    aux_2_2 = data_2_2_2.Fecha;
    for i = 1:length(names_ADME)
        index_2{i} = find(strcmp(aux_2_1,names_ADME{i})); % I save in a cell all the indices corresponding to names_ADME{i}.
        dates_2{i} = aux_2_2(index_2{i}); % I save in a cell all the dates corresponding to names_ADME{i}.
        aux = {};
        for j = 1:length(dates_2{i})
            aux{j} = datetime(dates_2{i}{j},'InputFormat','dd-MM-yyyy HH:mm:ss');
        end
        dates_time_2{i} = aux; % I transform the dates from strings to time variables.
        
        for j = 1:length(dates_time_2{i})
            if date_time >= dates_time_2{i}{j}
                instTurbine(i) = data_2_2_2.NU_Inst1(index_2{i}(j));
                MPTurbine(i) = data_2_2_2.NU_EnMP1(index_2{i}(j));
            end
        end
    end
    
    data_2(1) = CTR_Cost;
    data_2(2) = maxTurbine(2)*(instTurbine(2)-MPTurbine(2));
    data_2(3) = Mot_Cost;
    data_2(4) = maxTurbine(3)*(instTurbine(3)-MPTurbine(3));
    data_2(5) = PTA6_Cost;
    data_2(6) = maxTurbine(5)*(instTurbine(5)-MPTurbine(5));
    data_2(7) = PTB_Cost;
    data_2(8) = maxTurbine(1)*(instTurbine(1)-MPTurbine(1));
    data_2(9) = PTA78_Cost;
    data_2(10) = maxTurbine(4)*(instTurbine(4)-MPTurbine(4));

    %%%%%%%%%% Data #3:
    
    data_3_1 = xlsread(['./Data_Valores/Valores_de_Agua_ADME.xlsx'],1); % Bonete.
    data_3_2 = xlsread(['./Data_Valores/Valores_de_Agua_ADME.xlsx'],2); % Baygorria.
    data_3_3 = xlsread(['./Data_Valores/Valores_de_Agua_ADME.xlsx'],3); % Palmar.
    data_3_4 = xlsread(['./Data_Valores/Valores_de_Agua_ADME.xlsx'],4); % Salto Grande.
    
    data_3(1) = data_3_4(data_3_4==excelDate,2);
    data_3(2) = data_3_3(data_3_3==excelDate,2);
    data_3(3) = data_3_1(data_3_1==excelDate,2);
    data_3(4) = data_3_2(data_3_2==excelDate,2);
    
    %%%%%%%%%% Data #4:
    
    data_4_1 = readtable(['../Python/Represas_Data_2/Bonete_Data.xls']);
    
    row = find(contains(data_4_1.Var1,dateString_3));
    try
        data_4(1) = data_4_1.VolumenInicialEmbalsado(row); % Volume initial in hm^3.
        data_4(2) = data_4_1.AporteTe_rico(row); % Daily inflow of water in hm^3.
        data_4(3) = data_4_1.VolumenVertido(row); % Daily spilled water in hm^3.
    catch
        warning('Date not found for Bonete. See Data 4.');
    end
    
    %%%%%%%%%% Data #5:
    
    data_5_1 = readtable(['../Python/Represas_Data_2/Baygorria_Data.xls']);
    
    row = find(contains(data_5_1.Var1,dateString_3));
    try
        data_5(1) = data_5_1.VolumenInicialEmbalsado(row); % Volume initial in hm^3.
        data_5(2) = data_5_1.AporteTe_rico(row); % Daily inflow of water in hm^3.
        data_5(3) = data_5_1.VolumenVertido(row); % Daily spilled water in hm^3.
        
        row = find(contains(data_4_1.Var1,dateString_3_before));
        data_5(2) = data_5(2) - (6/24)*(data_4_1.VolumenTurbinado(row)+data_4_1.VolumenVertido(row));
        row = find(contains(data_4_1.Var1,dateString_3));
        data_5(2) = data_5(2) - (18/24)*(data_4_1.VolumenTurbinado(row)+data_4_1.VolumenVertido(row));
    catch
        warning('Date not found for Baygorria. See Data 5.');
    end
    
    %%%%%%%%%% Data #6:
    
    data_6_1 = readtable(['../Python/Represas_Data_2/Palmar_Data.xls']);
    
    row = find(contains(data_6_1.Var1,dateString_3));
    try
        data_6(1) = data_6_1.VolumenInicialEmbalsado(row); % Volume initial in hm^3.
        data_6(2) = data_6_1.AporteTe_rico(row); % Daily inflow of water in hm^3.
        data_6(3) = data_6_1.VolumenVertido(row); % Daily spilled water in hm^3.
        
        row = find(contains(data_5_1.Var1,dateString_3_before));
        data_6(2) = data_6(2) - (10/24)*(data_5_1.VolumenTurbinado(row)+data_5_1.VolumenVertido(row));
        row = find(contains(data_5_1.Var1,dateString_3));
        data_6(2) = data_6(2) - (14/24)*(data_5_1.VolumenTurbinado(row)+data_5_1.VolumenVertido(row));
    catch
        warning('Date not found for Palmar. See Data 6.');
    end
    
    %%%%%%%%%% Date #7:
    
    data_7_1 = readtable(['../Python/Represas_Data_2/Salto_Grande_Data_Horaria.csv']);
    data_7_1_1 = data_7_1.Var1; % Date.
    data_7_1_2 = data_7_1.Var2; % Hour.
    data_7_1_3 = data_7_1.CotaAguasArriba_m_;
    data_7_1_4 = data_7_1.CotaAguasAbajo_m_;
    data_7_1_5 = data_7_1.Salto_m_;
    
    row = find(data_7_1.Var1 == date_time);
    data_7(2) = data_7_1.CotaAguasArriba_m_(row(1)); % Initial level in m.
    
    data_7_2 = readtable(['../Python/Represas_Data_2/SG_Data/',dateString_4,'.csv']);
    % data_7_2 contains a lot of data with complicated characters.
    % Also, it is a table.
    aux = data_7_2.CaudalDeAporte{1}; % This is the string I want.
    % It has inflow from the river.
    aux_3 = strfind(aux,'.'); % I try to find the decimal point.
    % In this data, the decimal point indicates the first thousand.
    if isempty(aux_3)
        aux_2 =regexp(aux,'\d+(\.)?(\d+)?','match');
        data_7_2_1 = str2double([aux_2{:}]);
        % The last two lines extract the amount as a float number.
    else
        aux_2 =regexp(aux,'\d+(\.)?(\d+)?','match');
        data_7_2_1 = str2double([aux_2{:}]);
        % The last two lines extract the amount as a float number.
        data_7_2_1 = data_7_2_1*1000; % We correct the decimal point, and pass to m^3/s.
    end

    data_7(1) = data_7_2_1; % Average inflow in m^3/s.
    
    %%%%%%%%%% Date #8:
    
    data_8_1 = readtable(['../Python/Represas_Data_2/Bonete_Data_Horaria.csv']);
    data_8_1_1 = data_8_1.Var1; % Date.
    data_8_1_2 = data_8_1.Var2; % Hour.
    data_8_1_3 = data_8_1.CotaAguasArriba_m_;
    data_8_1_4 = data_8_1.CotaAguasAbajo_m_;
    data_8_1_5 = data_8_1.Salto_m_;
    
    row = find(data_8_1.Var1 == date_time);
    data_8_1_6 = data_8_1.CotaAguasAbajo_m_(row); % Level downstream of Bonete over time in m.
    data_8_1_6(25) = data_8_1_6(24); % I repeat the last value.
    
    data_8_2 = readtable(['../Python/Represas_Data_2/Baygorria_Data_Horaria.csv']);
    data_8_2_1 = data_8_2.Var1; % Date.
    data_8_2_2 = data_8_2.Var2; % Hour.
    data_8_2_3 = data_8_2.CotaAguasArriba_m_;
    data_8_2_4 = data_8_2.CotaAguasAbajo_m_;
    data_8_2_5 = data_8_2.Salto_m_;
    
    row = find(data_8_2.Var1 == date_time);
    data_8_2_6 = data_8_2.CotaAguasArriba_m_(row); % Level upstream of Baygorria in m.
    data_8_2_6(25) = data_8_2_6(24); % I repeat the last value.
    data_8_2_7 = data_8_2.CotaAguasAbajo_m_(row); % Level downstream of Baygorria in m.
    data_8_2_7(25) = data_8_2_7(24); % I repeat the last value.
    
    data_8_3 = readtable(['../Python/Represas_Data_2/Palmar_Data_Horaria.csv']);
    data_8_3_1 = data_8_3.Var1; % Date.
    data_8_3_2 = data_8_3.Var2; % Hour.
    data_8_3_3 = data_8_3.CotaAguasArriba_m_;
    data_8_3_4 = data_8_3.CotaAguasAbajo_m_;
    data_8_3_5 = data_8_3.Salto_m_;
    
    row = find(data_8_3.Var1 == date_time);
    data_8_3_6 = data_8_3.CotaAguasArriba_m_(row); % Level upstream of Palmar in m.
    data_8_3_6(25) = data_8_3_6(24); % I repeat the last value.
    data_8_3_7 = data_8_3.CotaAguasAbajo_m_(row); % Level downstream of Palmar in m.
    data_8_3_7(25) = data_8_3_7(24); % I repeat the last value.
    
    data_8_4 = readtable(['../Python/Represas_Data_2/Salto_Grande_Data_Horaria.csv']);
    data_8_4_1 = data_8_4.Var1; % Date.
    data_8_4_2 = data_8_4.Var2; % Hour.
    data_8_4_3 = data_8_4.CotaAguasArriba_m_;
    data_8_4_4 = data_8_4.CotaAguasAbajo_m_;
    data_8_4_5 = data_8_4.Salto_m_;
    
    row = find(data_8_4.Var1 == date_time);
    data_8_4_6 = data_8_4.CotaAguasArriba_m_(row); % Level upstream of SG in m.
    data_8_4_6(25) = data_8_4_6(24); % I repeat the last value.
    data_8_4_7 = data_8_4.CotaAguasAbajo_m_(row); % Level downstream of SG in m.
    data_8_4_7(25) = data_8_4_7(24); % I repeat the last value.
    
    data_8(:,1) = data_8_3_7;
    data_8(:,2) = data_8_4_7;
    data_8(:,7) = data_8_1_6-data_8_2_6;
    data_8(:,8) = data_8_2_7-data_8_3_6;
    
    %%%%%%%%%% Date #9:
    
    data_9_1 = readtable(['../Python/Represas_Data_2/Bonete_Data.xls']);
    data_9_2 = readtable(['../Python/Represas_Data_2/Baygorria_Data.xls']);
    
    row = find(contains(data_9_1.Var1,dateString_3_before));
    
    try
        aux = ones(61,1);
        data_9(:,1) = aux*(data_9_1.VolumenTurbinado(row)+data_9_1.VolumenVertido(row))*(1e6)/(24*3600);
        data_9(:,2) = aux*(data_9_2.VolumenTurbinado(row)+data_9_2.VolumenVertido(row))*(1e6)/(24*3600);
    catch
        warning('Date not found for Bonete or Baygorria. See Data 9.');
    end
    
    % What we did here, is to compute the averaged (turbine + spillage)
    % flow from the upstream dam, in a way to have the virtual control in
    % m^3/s for the next dams.
    
    %%%%%%%%%% Matrix:
    
    Matrix{1} = data_1;
    Matrix{2} = data_2;
    Matrix{3} = data_3;
    Matrix{4} = data_4;
    Matrix{5} = data_5;
    Matrix{6} = data_6;
    Matrix{7} = data_7;
    Matrix{8} = data_8;
    Matrix{9} = data_9;
    
    if 0 == exist('Historical','dir')
    	mkdir('Historical');
    end
    if 0 == exist('Historical/dailyData','dir')
    	mkdir('Historical/dailyData');
    end
    save([pwd '/Historical/dailyData/Day_',dateString_1,'.mat'],'Matrix');
    
end