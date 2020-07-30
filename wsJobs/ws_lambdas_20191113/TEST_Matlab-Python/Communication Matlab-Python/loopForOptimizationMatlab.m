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

while done == 0

    pause(2);
    try
        table = csvread('loopForOpt.csv');
    catch
        disp('loopForOpt.csv reading error!');
    end
    
    try
        if table(6) == 0
            x1 = table(1);
            x2 = table(2);
            x = [x1,x2];
            [f,g] = weirdCone(x);
            table(3) = f;
            table(4) = g(1);
            table(5) = g(2);
            table(6) = 1;
            csvwrite('loopForOpt.csv',table);
            disp('Work done!');
        else
            disp(['Table read, no FLAG = 0 found! Counter = ',num2str(counter)]);
        end
        counter = counter + 1;
    catch
        disp('Error at verifying condition!');
    end
    
end