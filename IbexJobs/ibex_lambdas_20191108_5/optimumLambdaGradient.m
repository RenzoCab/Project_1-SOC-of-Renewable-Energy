function output = optimumLambdaGradient(Day,optiNum)
    
    start = tic;
    
    figure(500);
    title('Duality GAP at each Iteration');
    xlabel('Number of Iterations');
    ylabel('Primal and Dual');
    grid minor; hold on;
    counter = 0;
    allTheData = {};
    
    % To test, I use T = 5 and S = 1.
    
    Exp_T = 6;
    Exp_S = 2;
    discLambda21 = 2;
    discLambda32 = 2;
    
    x0 = (2) * rand(1,discLambda21+discLambda32);
    Tol = 1e-4; % The value function units are USD 1e5.

    grad = [];
    hessian = [];
    
    if optiNum == 1
        options = optimoptions('fminsearch','Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol);
        [x,fval,exitflag,output] = fminsearch(@(x) Oracle(x),x0,options);
    elseif optiNum == 2
%     	options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'CheckGradients',true,'Display','notify',...
%             'PlotFcns',@optimplotfval,'TolFun',Tol,'FiniteDifferenceType','central','Algorithm','trust-region');
%         options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'CheckGradients',true,'Display','notify',...
%             'PlotFcns',@optimplotfval,'TolFun',Tol,'FiniteDifferenceType','forward','Algorithm','trust-region');
        options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol,'Algorithm','trust-region');
        [x,fval,exitflag,output,grad,hessian] = fminunc(@(x) Oracle(x),x0,options);
    elseif optiNum == 3
        options = optimoptions('patternsearch','Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol);
        [x,fval,exitflag,output] = patternsearch(@(x) Oracle(x),x0,options);
    elseif optiNum == 4
        options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol,'Algorithm','trust-region-reflective');
        [x,fval,exitflag,output] = fmincon(@(x) Oracle(x),x0,options);
    elseif optiNum == 5
        options = optimoptions('ga','Display','final',...
            'PlotFcns',@optimplotfval,'TolFun',Tol);
        [x,fval,exitflag,output] = ga(@(x) Oracle(x),length(x0),options);     
    end
    
    set(0,'CurrentFigure',500);
    saveas(gcf,[pwd '/Simulations/',num2str(optiNum),'_DG_',Day],'epsc');
    
    output = {x,fval,exitflag,output,grad,hessian,Exp_T,Exp_S,Day,Tol};
    save([pwd '/Simulations/',num2str(optiNum),'output_',Day,'.mat'],'output');
    
    for j = 1:discLambda21
        lambda21(j) = x(j);
    end
    for j = 1:discLambda32
        lambda32(j+discLambda21) = x(j+discLambda21);
    end
    
    save([pwd '/','Simulations/',num2str(optiNum),'allTheData_',Day,'.mat'],'allTheData');
    
    % After we have the optimal Lambda, we can use it in a more refined grid.
    Exp_S = 2;
    Exp_T = 8;
    Optimal_Solution(Day,3,1,0,...
        1,1,1,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
    
    toc(start)

    function [fx,gx] = Oracle(x)
        
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
    
        valDual = out(1);
        set(0,'CurrentFigure',500);
        plot(counter,valPrimal,'r*');
        plot(counter,valDual,'b*');
        plot(counter,valDualHJB,'k*');
        legend({'Primal OP','Dual OP','HJB dual problem'},'Location','eastoutside');
    
        fx = -valDualHJB; %-out(1);
        gx = -out(2:end);
    
        allTheData{counter,1} = valDualHJB;
        allTheData{counter,2} = valPrimal;
        allTheData{counter,3} = out;
        allTheData{counter,4} = x;
        disp(['Actual value: ',num2str(fx)]);
        
        if (optiNum == 1) || (optiNum == 3) || (optiNum == 5)
            gx = [];
        end
    
    end
    
end
