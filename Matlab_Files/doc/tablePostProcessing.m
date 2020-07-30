% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 09/12/2019

close all;
clear all;
clc;
set(0,'DefaultFigureVisible','off');

fSize = 12; % Font size for the plots.
Day = '20180302';
% finalDate = '20190101';
% Day = '20190101';
finalDate = '20200101';
count = 1;
date = {};
dateData = {};
expT = 6;
expS = 2;
normF = 1000;
costNorm = 1;
maxAbs = 1e5;
onlyClean = 0;
restartAll = 0;

while not(strcmp(Day,finalDate))
    
    if restartAll == 1
        delete([pwd '/Simulations/dataFrom_2_',Day,'.mat']);
        delete([pwd '/Simulations/dataFrom_4_',Day,'.mat']);
        delete([pwd '/Simulations/dataFrom_5_',Day,'.mat']);
        delete([pwd '/Simulations/dataFrom_6_',Day,'.mat']);
        delete([pwd '/Simulations/dataFrom_7_',Day,'.mat']);
    end
    
    disp('================================================================');
    disp(['=========================== ',Day,' ===========================']);
    disp('================================================================');
    
    disp('First part...');
    
    if 0 ~= exist([pwd '/Simulations/table_',Day,'.mat']) && 0 ~= exist([pwd '/Simulations/DG_',Day,'.eps'])
        
        load([pwd '/Simulations/table_',Day,'.mat']);
        
        if sum(size(saveTable) == [5 10]) == 2
        
            set(0,'DefaultFigureVisible','off');
            
            lambda21 = saveTable(5,1:2);
            lambda32 = saveTable(5,3:4);
            fx = -saveTable(5,5)*normF*costNorm; % Value function for the dual problem.
            % We multiply by normF because in loopForOptimizationMatlab and
            % oracleBlackBox we normalize using that value. Also, we change the
            % sign because before we were minimizing, and this is a concave
            % function.
            
            if onlyClean == 0
                
                if 0 ~= exist([pwd '/Simulations/dataFrom_2_',Day,'.mat'])
                    load([pwd '/Simulations/dataFrom_2_',Day,'.mat']);
                    valPrimal = dataFrom_2{1};
                else
                    close all;
                    set(0,'DefaultFigureVisible','off');
                    valPrimal = Optimal_Solution(Day,2,1,0,...
                        1,1,1,expS,expT,1,1,...
                        lambda21,lambda32,Day,1);
                end
               
                if 0 ~= exist([pwd '/Simulations/dataFrom_4_',Day,'.mat'])
                    load([pwd '/Simulations/dataFrom_4_',Day,'.mat']);
                    valSoft = dataFrom_4{1};
                else
                    close all;
                    set(0,'DefaultFigureVisible','off');    
                    valSoft = Optimal_Solution(Day,4,1,0,...
                        1,1,1,expS,expT,1,1,...
                        lambda21,lambda32,Day,1);
                end
                
                if 0 ~= exist([pwd '/Simulations/dataFrom_5_',Day,'.mat'])
                    load([pwd '/Simulations/dataFrom_5_',Day,'.mat']);
                    valBayFixed = dataFrom_5{1};
                else
                    close all;
                    set(0,'DefaultFigureVisible','off');    
                    valBayFixed = Optimal_Solution(Day,5,1,0,...
                        1,1,1,expS,expT,1,1,...
                        lambda21,lambda32,Day,1);
                end
                
                if 0 ~= exist([pwd '/Simulations/dataFrom_6_',Day,'.mat'])
                    load([pwd '/Simulations/dataFrom_6_',Day,'.mat']);
                    valDual = dataFrom_6{2};
                else
                    close all;
                    set(0,'DefaultFigureVisible','off');
                    valDual = Optimal_Solution(Day,6,1,0,...
                        1,1,1,expS,expT,1,1,...
                        lambda21,lambda32,Day,1);
                    valDual = valDual(1);
                end
                
                if 0 ~= exist([pwd '/Simulations/dataFrom_7_',Day,'.mat'])
                    load([pwd '/Simulations/dataFrom_7_',Day,'.mat']);
                    valReal = dataFrom_7{1};
                else
                    close all;
                    set(0,'DefaultFigureVisible','off');
                    valReal = Optimal_Solution(Day,7,1,0,...
                        1,1,1,expS,expT,1,1,...
                        lambda21,lambda32,Day,1);
                    valReal = valReal(1);
                end

                allFx(count)        = fx;
                allValPrimal(count) = valPrimal;
                allValDual(count)   = valDual;
                allBayFixed(count)  = valBayFixed;
                allValReal(count)   = valReal;
                allValSoft(count)   = valSoft;

                relativePrimalDG(count) = abs(valPrimal-fx)/abs(valPrimal);
                absolutePrimalDG(count) = abs(valPrimal-fx);
                relativeDualDG(count)   = abs(valPrimal-valDual)/abs(valPrimal);
                absoluteDualDG(count)   = abs(valPrimal-valDual);

                if relativePrimalDG(count) > 1
                    relativePrimalDG(count) = 1;
                end
                if absolutePrimalDG(count) > maxAbs
                    absolutePrimalDG(count) = maxAbs;
                end
                if relativeDualDG(count) > 1
                    relativeDualDG(count) = 1;
                end
                if absoluteDualDG(count) > maxAbs
                    absoluteDualDG(count) = maxAbs;
                end

                date{end + 1}     = Day;
                dateData{end + 1} = saveTable;

            end
        
        else
            
            relativePrimalDG(count) = inf;
            absolutePrimalDG(count) = inf;
            relativeDualDG(count) = inf;
            absoluteDualDG(count) = inf;
            delete([pwd '/Simulations/table_',Day,'.mat']);
            delete([pwd '/Simulations/DG_',Day,'.eps']);
            delete([pwd '/Simulations/table_2_',Day,'.mat']);
            delete([pwd '/Simulations/DG_2_',Day,'.eps']);
            
        end
                
    else
        
        relativePrimalDG(count) = inf;
        absolutePrimalDG(count) = inf;
        relativeDualDG(count) = inf;
        absoluteDualDG(count) = inf;
        delete([pwd '/Simulations/table_',Day,'.mat']);
        delete([pwd '/Simulations/DG_',Day,'.eps']);
        delete([pwd '/Simulations/table_2_',Day,'.mat']);
        delete([pwd '/Simulations/DG_2_',Day,'.eps']);
        
    end
    
    disp('Done!!!');
    disp('Second part...');
    
    if 0 ~= exist([pwd '/Simulations/table_2_',Day,'.mat']) && 0 ~= exist([pwd '/Simulations/DG_2_',Day,'.eps'])
            
        set(0,'DefaultFigureVisible','off');

        load([pwd '/Simulations/table_2_',Day,'.mat']);
        lambda21_2 = saveTable(5,1:4);
        lambda32_2 = saveTable(5,5:8);
        fx = -saveTable(5,9)*normF*costNorm; % Value function for the dual problem.
        % We multiply by normF because in loopForOptimizationMatlab and
        % oracleBlackBox we normalize using that value. Also, we change the
        % sign because before we were minimizing, and this is a concave
        % function.
        
        if onlyClean == 0
        
            if 0 == 1
            
                close all;
                set(0,'DefaultFigureVisible','off');    
                valSoft = Optimal_Solution(Day,4,1,0,...
                    0,0,0,expS,expT,1,1,...
                    lambda21_2,lambda32_2,Day,555);
                
                close all;
                set(0,'DefaultFigureVisible','off');    
                valSoft = Optimal_Solution(Day,5,1,0,...
                    0,0,0,expS,expT,1,1,...
                    lambda21_2,lambda32_2,Day,555);
                
                close all;
                set(0,'DefaultFigureVisible','off');
                valReal = Optimal_Solution(Day,7,1,0,...
                    0,0,0,expS,expT,1,1,...
                    lambda21_2,lambda32_2,Day,555);
                valReal = valReal(1);
                
            end
                
            close all;
            set(0,'DefaultFigureVisible','off');
            valPrimal = Optimal_Solution(Day,2,1,0,...
                0,0,0,expS,expT,1,1,...
                lambda21_2,lambda32_2,Day,555);
            
            close all;
            set(0,'DefaultFigureVisible','off');
            valDual = Optimal_Solution(Day,6,1,0,...
            	0,0,0,expS,expT,1,1,...
            	lambda21_2,lambda32_2,Day,555);
            valDual = valDual(1);
            
            allFx_2(count) = fx;
            allValPrimal_2(count) = valPrimal;
            allValDual_2(count) = valDual;
            allValReal_2(count) = valReal;
            allValSoft_2(count) = valSoft;

            relativePrimalDG_2(count) = abs(valPrimal-fx)/abs(valPrimal);
            absolutePrimalDG_2(count) = abs(valPrimal-fx);
            relativeDualDG_2(count) = abs(valPrimal-valDual)/abs(valPrimal);
            absoluteDualDG_2(count) = abs(valPrimal-valDual);

            if relativePrimalDG_2(count) > 1
                relativePrimalDG_2(count) = 1;
            end
            if absolutePrimalDG_2(count) > maxAbs
                absolutePrimalDG_2(count) = maxAbs;
            end
            if relativeDualDG_2(count) > 1
                relativeDualDG_2(count) = 1;
            end
            if absoluteDualDG_2(count) > maxAbs
                absoluteDualDG_2(count) = maxAbs;
            end

        end
        
    else
        
        delete([pwd '/Simulations/table_2_',Day,'.mat']);
        delete([pwd '/Simulations/DG_2_',Day,'.eps']);
        relativePrimalDG_2(count) = inf;
        absolutePrimalDG_2(count) = inf;
        relativeDualDG_2(count) = inf;
        absoluteDualDG_2(count) = inf;
        
    end
        
    disp('Done!!!');
    
	Day = nextDay(Day);
    count = count + 1;
    
end

if onlyClean == 0

    close all;
    set(0,'DefaultFigureVisible','on');
    pause(1);

    fin = length(dateData);
    days = 1:1:fin;

    figure;
    plot(days,relativePrimalDG(1:fin),'o');
    hold on;
    plot(days,relativePrimalDG_2(1:fin),'*');
    plot(days,0.1*ones(1,length(days)),'--r');
    title('Relative Duality GAP (OP-DHJB)');
    grid minor;
    xlabel('Days');
    xlim([days(1), days(end)]);
    ylim([0, 1]);
    h = legend('After 5 iterations with $\lambda_{21},\lambda_{32}\in R^2$','After 5 more iterations with $\lambda_{21},\lambda_{32}\in R^4$');
    set(h,'interpreter', 'latex');
    set(gca,'FontSize',fSize);
    saveas(gcf,[pwd '/Simulations/relativeDG'],'epsc');

    figure;
    plot(days,absolutePrimalDG(1:fin),'o');
    hold on;
    plot(days,absolutePrimalDG_2(1:fin),'*');
    title('Absolute Duality GAP (OP-DHJB)');
    grid minor;
    xlabel('Days');
    xlim([days(1), days(end)]);
    ylim([0, maxAbs]);
    ylabel('USD');
    h = legend('After 5 iterations with $\lambda_{21},\lambda_{32}\in R^2$','After 5 more iterations with $\lambda_{21},\lambda_{32}\in R^4$');
    set(h,'interpreter', 'latex');
    set(gca,'FontSize',fSize);
    saveas(gcf,[pwd '/Simulations/absoluteDG'],'epsc');

    figure;
    plot(days,relativeDualDG(1:fin),'o');
    hold on;
    plot(days,relativeDualDG_2(1:fin),'*');
    plot(days,0.1*ones(1,length(days)),'--r');
    title('Relative Duality GAP (OP-DOP)');
    grid minor;
    xlabel('Days');
    xlim([days(1), days(end)]);
    ylim([0, 1]);
    h = legend('After 5 iterations with $\lambda_{21},\lambda_{32}\in R^2$','After 5 more iterations with $\lambda_{21},\lambda_{32}\in R^4$');
    set(h,'interpreter', 'latex');
    set(gca,'FontSize',fSize);
    saveas(gcf,[pwd '/Simulations/relativeDualDG'],'epsc');

    figure;
    plot(days,absoluteDualDG(1:fin),'o');
    hold on;
    plot(days,absoluteDualDG_2(1:fin),'*');
    title('Absolute Duality GAP (OP-DOP)');
    grid minor;
    xlabel('Days');
    xlim([days(1), days(end)]);
    ylim([0, maxAbs]);
    ylabel('USD');
    h = legend('After 5 iterations with $\lambda_{21},\lambda_{32}\in R^2$','After 5 more iterations with $\lambda_{21},\lambda_{32}\in R^4$');
    set(h,'interpreter', 'latex');
    set(gca,'FontSize',fSize);
    saveas(gcf,[pwd '/Simulations/absoluteDualDG'],'epsc');

    figure;
    % allFx(count) = fx;
    % allValPrimal(count) = valPrimal;
    % allValDual(count) = valDual;
    % allValReal(count) = valReal;
    plot(days,allValReal(1:fin),'*');
    hold on;
    plot(days,allValSoft(1:fin),'o','MarkerSize',10);
    plot(days,allValPrimal(1:fin),'o','MarkerSize',10);
    plot(days,allValDual(1:fin),'o');
    plot(days,allFx(1:fin),'*');
    title('All Costs Comparison');
    grid minor;
    xlabel('Days');
    xlim([days(1), days(end)]);
    ylabel('USD');
    legend('Historical Cost','Soft Primal Solution','Admissible Primal Solution (OP)','Dual Cost-to-Go (DOP)','Dual HJB (DHJB)');
    set(gca,'FontSize',fSize);
    saveas(gcf,[pwd '/Simulations/allCostsComparison'],'epsc');

    accumReal       = allValReal(1);
    accumPrimal     = allValPrimal(1);
    accumSoft       = allValSoft(1);
    accumBayFixed   = allBayFixed(1);
    for i = 2:fin
        if allValReal(i) ~= Inf
            accumReal(i)     = accumReal(i-1)     + allValReal(i);
            accumPrimal(i)   = accumPrimal(i-1)   + allValPrimal(i);
            accumSoft(i)     = accumSoft(i-1)     + allValSoft(i);
            accumBayFixed(i) = accumBayFixed(i-1) + allBayFixed(i);
        end
    end
    figure;
    plot(days,accumReal,'*--');
    hold on;
    plot(days,accumSoft,'o--','MarkerSize',10);
    plot(days,accumPrimal,'*--');
    plot(days,accumBayFixed,'*--');
    title('Accumulated Costs');
    grid minor;
    xlabel('Days');
    xlim([days(1), days(end)]);
    ylabel('USD');
    legend('Real Accumulated Cost','Soft (Addmisible) Accumulated Cost','Primal (Addmisible) Accumulated Cost','No Se');
    set(gca,'FontSize',fSize);
    saveas(gcf,[pwd '/Simulations/allAccumulatedCosts'],'epsc');

end

%% Print Tables:

if onlyClean == 0

    for i = 1:fin
        if relativePrimalDG(i) > 0.2 && relativePrimalDG(i) ~= Inf

            disp('==================== More 20 % ====================');
            disp(date{i});
            disp(dateData{i});
            disp('==================== More 20 % ====================');

        elseif relativePrimalDG(i) > 0.1 && relativePrimalDG(i) ~= Inf

            disp('==================== More 10 % ====================');
            disp(date{i});
            disp(dateData{i});
            disp('==================== More 10 % ====================');

        end
    end

end