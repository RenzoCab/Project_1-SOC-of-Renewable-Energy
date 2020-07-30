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

Day = '20180302';
finalDate = '20180501';
count = 1;
expT = 6;
expS = 2;
normF = 1000;

while not(strcmp(Day,finalDate))
    
    if 0 ~= exist([pwd '/Simulations/table_',Day,'.mat'])
       
        load([pwd '/Simulations/table_',Day,'.mat']);
        lambda21 = saveTable(5,1:2);
        lambda32 = saveTable(5,3:4);
        fx = -saveTable(5,5)*normF; % Value function for the dual problem.
        % We multiply by normF because in loopForOptimizationMatlab and
        % oracleBlackBox we normalize using that value. Also, we change the
        % sign because before we were minimizing, and this is a concave
        % function.
                
        valPrimal = Optimal_Solution(Day,2,1,0,...
            0,0,0,expS,expT,1,1,...
            lambda21,lambda32,Day,1);
        
        allFx(count) = fx;
        allValPrimal(count) = valPrimal;
        relativeDG(count) = abs(valPrimal-fx)/abs(valPrimal);
        absoluteDG(count) = abs(valPrimal-fx);
                
    end
    
    if 0 ~= exist([pwd '/Simulations/table_2_',Day,'.mat'])
       
        load([pwd '/Simulations/table_2_',Day,'.mat']);
        lambda21_2 = saveTable(5,1:4);
        lambda32_2 = saveTable(5,5:8);
        fx = -saveTable(5,9)*normF; % Value function for the dual problem.
        % We multiply by normF because in loopForOptimizationMatlab and
        % oracleBlackBox we normalize using that value. Also, we change the
        % sign because before we were minimizing, and this is a concave
        % function.
                
        valPrimal = Optimal_Solution(Day,2,1,0,...
            0,0,0,expS,expT,1,1,...
            lambda21,lambda32,Day,1);
        
        allFx_2(count) = fx;
        allValPrimal_2(count) = valPrimal;
        relativeDG_2(count) = abs(valPrimal-fx)/abs(valPrimal);
        absoluteDG_2(count) = abs(valPrimal-fx);
        
        count = count + 1;
                
    elseif 0 ~= exist([pwd '/Simulations/table_',Day,'.mat'])
        
        relativeDG_2(count) = 0;
        absoluteDG_2(count) = 0;
        
        count = count + 1;
        
    end
        
        Day = nextDay(Day);
    
end

days = 1:1:length(relativeDG);
plot(days,relativeDG,'*-');
hold on;
plot(days,relativeDG_2,'*-');
title('Relative Duality GAP');
grid on;
xlabel('Days');
xlim([days(1), days(end)]);
saveas(gcf,[pwd '/Simulations/relativeDG'],'epsc');

figure;
plot(days,absoluteDG,'*-');
hold on;
plot(days,absoluteDG_2,'*-');
title('Absolute Duality GAP');
grid on;
xlabel('Days');
xlim([days(1), days(end)]);
saveas(gcf,[pwd '/Simulations/absoluteDG'],'epsc');
