% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 12/12/2019

Day = '20180302';
finalDate = '20190101';
% finalDate = '20180303';
normF = 1000;
costNorm = 1;
expT = 6;
expS = 2;

while not(strcmp(Day,finalDate))
   
    if 0 ~= exist([pwd '/Simulations/table_',Day,'.mat'])
        
        load([pwd '/Simulations/table_',Day,'.mat']);
        
        if sum(size(saveTable) == [5 10]) == 2
                    
            lambda21 = saveTable(5,1:2);
            lambda32 = saveTable(5,3:4);
            fx = -saveTable(5,5)*normF*costNorm; % Value function for the dual problem.
            % We multiply by normF because in loopForOptimizationMatlab and
            % oracleBlackBox we normalize using that value. Also, we change the
            % sign because before we were minimizing, and this is a concave
            % function.

%             close all;
%             set(0,'DefaultFigureVisible','off');
%             Optimal_Solution(Day,2,1,0,...
%                 1,1,1,expS,expT,1,1,...
%                 lambda21,lambda32,Day,1);

            close all;
            set(0,'DefaultFigureVisible','off');    
            Optimal_Solution(Day,4,1,0,...
                1,1,1,expS,expT,1,1,...
                lambda21,lambda32,Day,555);                

%             close all;
%             set(0,'DefaultFigureVisible','off');
%             Optimal_Solution(Day,6,1,0,...
%                 1,1,1,expS,expT,1,1,...
%                 lambda21,lambda32,Day,1);

            close all;
            set(0,'DefaultFigureVisible','off');
            Optimal_Solution(Day,7,1,0,...
                1,1,1,expS,expT,1,1,...
                lambda21,lambda32,Day,555);
                    
        end
    
    end
    
    Day = nextDay(Day);

    
end