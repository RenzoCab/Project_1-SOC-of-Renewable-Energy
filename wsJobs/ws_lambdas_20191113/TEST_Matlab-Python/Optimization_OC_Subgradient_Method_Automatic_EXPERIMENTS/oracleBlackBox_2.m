% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 05/12/2019

function output = oracleBlackBox_2(x0,Day,Exp_T,Exp_S,discLambda21,discLambda32,...
    counter,savePlotsIterations,normF,normG)
    
    disp('Oracle Starting...');
    start = tic;
    
    % We create Figure 500.
    figure(500);
    title('Duality GAP at each Iteration');
    xlabel('Number of Iterations');
    ylabel('Primal and Dual');
    grid minor; hold on;
    
    % Evaluation. We obtain -[dual function] and -[subgradient].
    output = Oracle(x0);
    
    % We save Figure 500.
    set(0,'CurrentFigure',500);
    saveas(gcf,[pwd '/Simulations/DG_2_',Day],'epsc');
        
    toc(start)

    function output_2 = Oracle(x)
                
        % We create the two Lambdas:
        for i = 1:discLambda21
            lambda21(i) = x(i);
        end
        for i = 1:discLambda32
            lambda32(i) = x(i+discLambda21);
        end
        
        % We evaluate the dual function using HJB:
        valDualHJB = Optimal_Solution(Day,1,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
        
        % I find the admissible optimal path, but with no penalization H^1.
        % This is an upper bound for the real primal solution:
        valPrimal = Optimal_Solution(Day,2,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
        
        % We find the dual optimal path. It is necessary to compute the subgradient:
        out = Optimal_Solution(Day,6,1,0,...
        savePlotsIterations,savePlotsIterations,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1); % gx is in m^3/s.
        if savePlotsIterations
            if 0 == exist([pwd '/Simulations/Simulation_2_',Day,'_',num2str(counter)],'dir')
                mkdir([pwd '/Simulations/Simulation_2_',Day,'_',num2str(counter)]);
            end
            for i = 101:133
                set(0,'CurrentFigure',i);
                saveas(gcf,[pwd '/Simulations/Simulation_2_',Day,'_',num2str(counter),'/',num2str(i)],'epsc');
            end
            for i = 1000:1004
                set(0,'CurrentFigure',i);
                saveas(gcf,[pwd '/Simulations/Simulation_2_',Day,'_',num2str(counter),'/',num2str(i)],'epsc');
            end
            i = 11;
            set(0,'CurrentFigure',i);
            saveas(gcf,[pwd '/Simulations/Simulation_2_',Day,'_',num2str(counter),'/',num2str(i)],'epsc');
        end
    
        % To plot the duality GAP:
        valDual = out(1);
        set(0,'CurrentFigure',500);
        plot(counter,valPrimal,'r*');
        plot(counter,valDual,'b*');
        plot(counter,valDualHJB,'k*');
        ylabel('USD');
        xlabel('Evaluation');
        legend({'Primal OP','Dual OP','HJB dual problem'},'Location','eastoutside');
    
        fx = -valDualHJB/normF; % -out(1);
        gx = -out(2:end)/normG;
        output_2 = [fx,gx];
        % We divide between normF and normG to make it friendlier for the optimization.
    
    end
    
end