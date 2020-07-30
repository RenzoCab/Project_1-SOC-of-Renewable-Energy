function output = optimumLambdaGradient(Day,optiNum)
    
    start = tic;
    
    testing = 1;
    if testing == 1
        vector = [-10,-5,0,5,10];
        startingPoint = {};
        for k = 1:5
            for j = 1:5
                startingPoint{end+1} = [vector(k)*ones(1,2),vector(j)*ones(1,2)];
            end
        end
    end
    
    figure(500);
    title('Duality GAP at each Iteration');
    xlabel('Number of Iterations');
    ylabel('Primal and Dual');
    grid minor; hold on;
    counter = 0;
    allTheData = {};
    numOfIterations = 0; % Artificial number of evaluations counter.
    limitOfIterations = 100; % Artificial maximum number of evaluations.
    maxFuncEval = 10; % Maximum number of function evaluations.
    maxIter = 5; % Maximum number of iterations.
    numberOfDivisions = 1; % Times we divide the Lagrange multiplier.
    
    Exp_T = 6;
    Exp_S = 2;
    discLambda21 = 2;
    discLambda32 = 2;
    
    tools = {'fminsearch','fminunc','patternsearch','fmincon','ga'};
    disp(['Optimization tool: ',tools{optiNum}]);
    
    x0 = (2) * rand(1,discLambda21+discLambda32); % <<< Rand or Ones?
    x0 = (2) * ones(1,discLambda21+discLambda32);
    Tol = 1e-4; % The value function units are USD 1e5.

    grad = [];
    hessian = [];
    
    if optiNum == 1
        
        options = optimset('Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol);
        [x,fval,exitflag,output] = fminsearch(@(x) Oracle(x),x0,options);
        
    elseif optiNum == 2
        
%     	options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'CheckGradients',true,'Display','notify',...
%             'PlotFcns',@optimplotfval,'TolFun',Tol,'FiniteDifferenceType','central','Algorithm','trust-region');
%         options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'CheckGradients',true,'Display','notify',...
%             'PlotFcns',@optimplotfval,'TolFun',Tol,'FiniteDifferenceType','forward','Algorithm','trust-region');
        options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol,'Algorithm','trust-region','OutputFcn', @outfun,...
            'MaxFunEvals',maxFuncEval,'MaxIterations',maxIter);
        swtch = 0; % With 0 or 1 we divide by turns, with 2 we divide both at once.
        for z = 1:numberOfDivisions
            if testing == 1
                for h = 1:length(startingPoint)
                    x0 = startingPoint{h};
                    [x,fval,exitflag,output,grad,hessian] = fminunc(@(x) Oracle(x),x0,options);
                    set(0,'CurrentFigure',500); clf(500); hold on;
                end
            else
                [x,fval,exitflag,output,grad,hessian] = fminunc(@(x) Oracle(x),x0,options);
            end
            xAux1 = x(1:length(x)/2);
            xAux2 = x(length(x)/2+1:end);
            if swtch == 0
                discLambda21 = discLambda21 * 2;
                xAux1 = repelem(xAux1,2);
                x0 = [xAux1,xAux2];
                swtch = 1;
            elseif swtch == 1
                discLambda32 = discLambda32 * 2;
                xAux2 = repelem(xAux2,2);
                x0 = [xAux1,xAux2];
                swtch = 0;
            elseif swtch == 2
                discLambda21 = discLambda21 * 2;
                discLambda32 = discLambda32 * 2;
                x0 = repelem(x,2);
            end
        end
        
    elseif optiNum == 3
        
        options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
        A       = [];
        b       = [];
        Aeq     = [];
        beq     = [];
        lb      = [];
        ub      = [];
        nonlcon = [];
        [x,fval,exitflag,output] = patternsearch(@(x) Oracle(x),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        
    elseif optiNum == 4
        
        options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','notify',...
            'PlotFcns',@optimplotfval,'TolFun',Tol,'Algorithm','trust-region-reflective');
        A       = [];
        b       = [];
        Aeq     = [];
        beq     = [];
        lb      = [];
        ub      = [];
        nonlcon = [];
        [x,fval,exitflag,output] = fmincon(@(x) Oracle(x),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        
    elseif optiNum == 5
        
        options = optimoptions('ga','Display','final',...
            'PlotFcns',@optimplotfval,'TolFun',Tol);
        [x,fval,exitflag,output] = ga(@(x) Oracle(x),length(x0),options);
    end
    
    set(0,'CurrentFigure',500);
    saveas(gcf,[pwd '/Simulations/',num2str(optiNum),'_DG_',Day],'epsc');
    
    output = {x,fval,exitflag,output,grad,hessian,Exp_T,Exp_S,Day,Tol};
    save([pwd '/Simulations/',num2str(optiNum),'_output_',Day,'.mat'],'output');
    
    for j = 1:discLambda21
        lambda21(j) = x(j);
    end
    for j = 1:discLambda32
        lambda32(j+discLambda21) = x(j+discLambda21);
    end
    
    save([pwd '/','Simulations/',num2str(optiNum),'_allTheData_',Day,'.mat'],'allTheData');
    
    % After we have the optimal Lambda, we can use it in a more refined grid.
    Exp_S = 3;
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
        ylabel('USD');
        xlabel('Evaluation');
        legend({'Primal OP','Dual OP','HJB dual problem'},'Location','eastoutside');
        if testing == 1
            saveas(gcf,[pwd '/Simulations/',num2str(h),'_DG_',Day],'epsc');
        else
            saveas(gcf,[pwd '/Simulations/',num2str(optiNum),'_DG_',Day],'epsc');
        end
    
        fx = -valDualHJB/1000; %-out(1);
        gx = -out(2:end)/10;
        % We divide between 10 and 1000 to make it friendlier for the optimization.
        
        allTheData{counter,1} = valDualHJB;
        allTheData{counter,2} = valPrimal;
        allTheData{counter,3} = out;
        allTheData{counter,4} = x;
        disp(['Actual value: ',num2str(fx)]);
        
        if (optiNum == 1) || (optiNum == 3) || (optiNum == 5)
            gx = [];
        end
    
    end

    function stop = outfun(x, optimValues, state)
        numOfIterations = numOfIterations + 1;
        disp(['numOfIterations = ',num2str(numOfIterations)]);
        if numOfIterations == limitOfIterations
            stop = true;
            cprintf('err','The optimization will STOP!!!\n');
            numOfIterations = 0;
        else
            stop = false;
        end
    end
    
end