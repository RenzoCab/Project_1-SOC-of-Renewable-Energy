function output = finiteDifferences(Day)
    
    start = tic;
    
    counter = 0;
    
    % To test, I use T = 5 and S = 1.
    
    Exp_T = 6;
    Exp_S = 2;
    discLambda21 = 2;
    discLambda32 = 2;
    
    x0 = 1 - 2*rand(1,discLambda21+discLambda32);
    
    delta = 0.1;
    
    allData = {};
    
    allData{1} = Oracle(x0 + [0 0 0 0]);
    allData{2} = Oracle(x0 + [0 0 0 delta]);
    allData{3} = Oracle(x0 + [0 0 0 -delta]);
    allData{4} = Oracle(x0 + [0 0 delta 0]);
    allData{5} = Oracle(x0 + [0 0 -delta 0]);
    allData{6} = Oracle(x0 + [0 delta 0 0]);
    allData{7} = Oracle(x0 + [0 -delta 0 0]);
    allData{8} = Oracle(x0 + [delta 0 0 0]);
    allData{9} = Oracle(x0 + [-delta 0 0 0]);
    allData{10} = delta;
    
    save([pwd '/','Simulations/FD_',Day,'.mat'],'allData');
    
    toc(start)

    function [output] = Oracle(x)
        
        counter = counter + 1;
        
        % We create the two Lambdas:
        for i = 1:discLambda21
            lambda21(i) = x(i);
        end
        for i = 1:discLambda32
            lambda32(i) = x(i+discLambda21);
        end
        
        savePlotsIterations = 1;

        valDualHJB = Optimal_Solution(Day,1,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
        
        % I find the admissible optimal path but with no penalization H^1.
        valPrimal = Optimal_Solution(Day,2,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
        
        out = Optimal_Solution(Day,6,1,0,...
        savePlotsIterations,savePlotsIterations,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1); % gx is in m^3/s.
        if savePlotsIterations
            if 0 == exist([pwd '/Simulations/Simulation_',Day,'_',num2str(counter)],'dir')
                mkdir([pwd '/Simulations/Simulation_',Day,'_',num2str(counter)]);
            end
            for i = 101:133
                set(0,'CurrentFigure',i);
                saveas(gcf,[pwd '/Simulations/Simulation_',Day,'_',num2str(counter),'/',num2str(i)],'epsc');
            end
            for i = 1000:1004
                set(0,'CurrentFigure',i);
                saveas(gcf,[pwd '/Simulations/Simulation_',Day,'_',num2str(counter),'/',num2str(i)],'epsc');
            end
            i = 11;
            set(0,'CurrentFigure',i);
            saveas(gcf,[pwd '/Simulations/Simulation_',Day,'_',num2str(counter),'/',num2str(i)],'epsc');
        end
        close all;
    
        valDual = out(1);
    
        fx = -valDualHJB; %-out(1);
        gx = -out(2:end);
    
        disp(['Actual value: ',num2str(fx)]);
        
        output = [fx,gx];
    
    end
    
end
