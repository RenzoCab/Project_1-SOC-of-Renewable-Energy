function [A] = convexityDualFunction(Day)
    
    start = tic;
    
    i = 0; % i is in {0,1,...,Exp_T-2}. We divide the support of the Lagrange Multiplayer into (3*2^i) partitions.
    Exp_T = 5;
    Exp_S = 1;
    k = 100;
    
    for j = -k:1:k
        smallVal = 1e-5;
        x0 = j*ones(1,3*2^i)*smallVal;
        A(j+k+1) = oracle(x0);
        B(j+k+1) = j*smallVal;
    end
    
    figure; plot(B,A);
    
    toc(start)
    
    function [fx] = oracle(x)

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

        fx = Optimal_Solution('Opt',1,1,0,...
        0,0,0,Exp_S,Exp_T,1,1,...
        Vec,0,Day,i);
    
%         ffx = Optimal_Solution('Opt',2,1,0,...
%         1,1,1,Exp_S,Exp_T,1,1,...
%         Vec,0,Day,i);

%         gx = Optimal_Solution('Opt',6,1,0,...
%         0,0,0,Exp_S,Exp_T,1,1,...
%         Vec,0,Day,i);
        
        disp(['Actual value: ',num2str(fx)]);
    
    end
    
end
