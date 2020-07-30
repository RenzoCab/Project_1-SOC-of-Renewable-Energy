function output = optimumLambdaGradient(Day)
    
    start = tic;
    
    figure(500);
    title('Duality GAP at each Iteration');
    xlabel('Number of Iterations');
    ylabel('Primal and Dual');
    grid minor; hold on;
    counter = 0;
    
    % To test, I use T = 5 and S = 1.
    
    Exp_T = 6;
    Exp_S = 2;
    discLambda21 = 2;
    discLambda32 = 2;
    
    x0 = (1) * ones(1,discLambda21+discLambda32);
    Tol = 1e-5; % The value function units are USD 1e5.

% 	options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'CheckGradients',true,'Display','notify',...
%         'PlotFcns',@optimplotfval,'TolFun',Tol,'FiniteDifferenceType','central','Algorithm','trust-region');
    options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','notify',...
        'PlotFcns',@optimplotfval,'TolFun',Tol,'Algorithm','trust-region');
    
    [x,fval,exitflag,output,grad,hessian] = fminunc(@(x) Oracle(x),x0,options);
    
    output = {x,fval,exitflag,output,grad,hessian,Exp_T,Exp_S,Day,Tol};
    save(['output_',Day,'.mat'],'output');
    
    for j = 1:discLambda21
        lambda21(j) = x(j);
    end
    for j = 1:discLambda32
        lambda32(j+discLambda21) = x(j+discLambda21);
    end
    
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

        valDualHJB = Optimal_Solution(Day,1,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
    
        % I find the admissible optimal path but with no penalization H^1.
        valPrimal = Optimal_Solution(Day,2,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1);
        
        out = Optimal_Solution(Day,6,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        lambda21,lambda32,Day,1); % gx is in m^3/s.
    
        valDual = out(1);
        set(0,'CurrentFigure',500);
        plot(counter,valPrimal,'r*');
        plot(counter,valDual,'b*');
        plot(counter,valDualHJB,'k*');
    
        fx = -valDualHJB; %-out(1);
        gx = -out(2:end);
    
        disp(['Actual value: ',num2str(fx)]);
    
    end
    
end
