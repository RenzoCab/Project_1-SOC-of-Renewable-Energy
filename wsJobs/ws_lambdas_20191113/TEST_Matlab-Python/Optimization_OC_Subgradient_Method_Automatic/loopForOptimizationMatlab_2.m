% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 05/12/2019

close all;
clear all;
clc;

parpool(feature('numcores'));
done = 0;
counter = 0;
counter_2 = 1;
tableCounter = 0;
initialDate = '00000000';
date = initialDate;
expT = 6;
expS = 2;
discL21 = 4;
discL32 = 4;
savePlots = 1;
normF = 1000;
normG = 10;

while done == 0

    pause(2);
    try
        table = csvread('loopForOpt_2.csv');
        dateLoad = csvread('date_2.csv');
        dateLoad = num2str(dateLoad);
        if not(strcmp(date,dateLoad))
            date = dateLoad;
            tableCounter = 0;
            close all;
            try
                auxDate = previousDay(date);
                save([pwd '/','Simulations/table_2_',auxDate,'.mat'],'saveTable');
            end
        end
    catch
        disp('loopForOpt_2.csv or date_2.csv reading error!');
    end
    
    % table = [x1,x2,x3,x4,x5,x6,x7,x8,fx,gx1,gx2,gx3,gx4,gx5,gx6,gx7,gx8,FLAG];
    
%     try
        if (table(end) == 0) && not(strcmp(date,initialDate))
            
            tableCounter = tableCounter + 1;
            
            x = table(1:8);
            
            output = oracleBlackBox_2(x,date,expT,expS,discL21,discL32,tableCounter,...
                savePlots,normF,normG);
            % output = [fx,gx], fx \in \R, gx \in \R^4.
            
            table(9:17) = output(1:9);
            table(end) = 1;
            csvwrite('loopForOpt_2.csv',table);
            disp('Work done!');
            
            saveTable(tableCounter,:) = table;
            
        else
            disp(['Table read, no FLAG = 0 found! Counter = ',num2str(counter)]);
        end
        counter = counter + 1;
%     catch
%         disp('Error at verifying condition!');
%     end
    
end