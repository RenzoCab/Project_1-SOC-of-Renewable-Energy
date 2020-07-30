% Author: Renzo Caballero
% Institution: KAUST (King Abdullah University of Science and Technology)
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 20/11/2019

close all;
clear all;
clc;

done = 0;
counter = 0;

while done == 0

    pause(2);
    try
        table = csvread('loop.csv');
    catch
        disp('loop.csv not found!');
    end
    
    try
        if table(end) == 0
            if rand < 0.2
                table(end) = 3;
            else
                table(end) = 1;
            end
            csvwrite('loop.csv',table);
            disp('Work done!');
        elseif table(end) == 2
            table(end) = 1;
            csvwrite('loop.csv',table);
            disp('Wow! Fortran sent a 2!');
        else
            disp(['Table read, no FLAG = 0 found! Counter = ',num2str(counter)]);
        end
        counter = counter + 1;
    catch
        disp('Error at verifying condition!');
    end
    
end