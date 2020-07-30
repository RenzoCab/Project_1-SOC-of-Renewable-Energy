% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 14/11/2019

close all;
clear all;
clc;

done = 0;
counter = 0;
tableCounter = 0;
day = '20180302';
expT = 6;
expS = 2;
discL21 = 2;
discL32 = 2;
savePlots = 1;
normF = 1000;
normG = 10;

while done == 0

    pause(2);
    try
        table = csvread('loopForOpt.csv');
    catch
        disp('loopForOpt.csv reading error!');
    end
    
    % table = [x1,x2,x3,x4,fx,gx1,gx2,gx3,gx4,FLAG];
    
%     try
        if table(10) == 0
            
            x1 = table(1);
            x2 = table(2);
            x3 = table(3);
            x4 = table(4);
            x = [x1,x2,x3,x4];
            
            output = oracleBlackBox(x,day,expT,expS,discL21,discL32,counter,...
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
            
            tableCounter = tableCounter + 1;
            saveTable(tableCounter,:) = table;
        else
            disp(['Table read, no FLAG = 0 found! Counter = ',num2str(counter)]);
        end
        counter = counter + 1;
%     catch
%         disp('Error at verifying condition!');
%     end
    
end