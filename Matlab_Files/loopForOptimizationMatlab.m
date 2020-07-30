% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 02/07/2020

close all;
clear all;
clc;

try
    parpool(feature('numcores'));
catch
    disp('Parpool already started...');
end

done = 0;
counter = 0;
counter_2 = 1;
tableCounter = 0;
initialDate = '00000000';
date = initialDate;
expT = 6;
expS = 2;
discL21 = 2;
discL32 = 2;
savePlots = 0;
normF = 1000;
normG = 10;

% date = '20190102';
% writematrix(date,'date.csv');
% load('date.csv');

while done == 0

    pause(2);
    try
        table = csvread('loopForOpt.csv');
        dateLoad = csvread('date.csv');
        dateLoad = num2str(dateLoad);
        if not(strcmp(date,dateLoad))
            date = dateLoad;
            tableCounter = 0;
            close all;
            try
                auxDate = previousDay(date);
                save([pwd '/','Simulations/table_',auxDate,'.mat'],'saveTable');
            end
        end
    catch
        disp('loopForOpt.csv or date.csv reading error!');
    end
    
    % table = [x1,x2,x3,x4,fx,gx1,gx2,gx3,gx4,FLAG];
    
%     try
        if (table(10) == 0) && not(strcmp(date,initialDate))
            
            tableCounter = tableCounter + 1;
            
            x1 = table(1);
            x2 = table(2);
            x3 = table(3);
            x4 = table(4);
            x = [x1,x2,x3,x4];
            
            output = oracleBlackBox(x,date,expT,expS,discL21,discL32,tableCounter,...
                savePlots,normF,normG);
            % output = [fx,gx], fx \in \R, gx \in \R^4.
            
            table(5) = output(1);
            table(6) = output(2);
            table(7) = output(3);
            table(8) = output(4);
            table(9) = output(5);
            table(10) = 1;
            csvwrite('loopForOpt.csv',table);
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