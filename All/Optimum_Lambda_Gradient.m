function [x,fval] = Optimum_Lambda_Gradient(Day)
    
    start = tic;
    
    i = 3; % i is in {0,1,...,Exp_T-2}
    Exp_T = 8;
    Exp_S = 2;
    x0 = 0*ones(1,3*2^i);
    Tol = 1e-9;

	options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','notify',...
        'PlotFcns',@optimplotfval,'TolFun',Tol); % Default: 'TolFun',1e-6.
    [x,fval,exitflag,output] = fminunc(@(x) Oracle_2(x),x0,options);
    output = {x,fval,exitflag,output,i,Exp_T,Exp_S,Day,Tol};
    h1 = gcf;
    h2 = figure(1003); % Change this variable for each figure.
    copyobj(get(h1,'children'),h2);
    grid minor;
    for m = 1000:1003
    	set(0,'CurrentFigure',m);
        if 0 == exist(['Simulation_Day_NB_',num2str(Day)],'dir')
            mkdir(['Simulation_Day_NB_',num2str(Day)]);
        end
     	saveas(gcf,[pwd '/',['Simulation_Day_NB_',num2str(Day)],'/Opt_',num2str(m)],'epsc');
    end
    save([pwd '/Simulation_Day_NB_',num2str(Day),'/Lambda.mat'],'output');
    
    toc(start)
    
    function [fx,gx] = Oracle_2(x)

        if exist('Simulation_Opt.mat','file') == 2
            delete Simulation_Opt.mat
        end

        T = 2^Exp_T;
        Z21 = 2^Exp_T/4;
        Delta = (T-Z21)/(3*2^i);
        for n = 0:length(x)-1
            aux_x(1+n*Delta:(n+1)*Delta) = x(n+1);
        end
        aux_x(end+1) = aux_x(end);
        Vec = [zeros(1,Z21),aux_x,zeros(1,Z21)];

%         fx = -1*Admissible_Solution_6('Opt',1,1,0,...
%         0,0,0,Exp_S,Exp_T,1,1,...
%         Vec,0,Day,i);
% 
%         gx = -1*Admissible_Solution_6('Opt',6,1,0,...
%         0,0,0,Exp_S,Exp_T,1,1,...
%         Vec,0,Day,i);
        
        fx = -1*Admissible_Solution_7('Opt',1,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        Vec,0,Day,i);

        gx = -1*Admissible_Solution_7('Opt',6,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        Vec,0,Day,i);
    
    end
    
end