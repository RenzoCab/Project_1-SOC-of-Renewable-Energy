function [fx,gx] = Oracle(x)

    if exist('Simulation_Opt.mat','file') == 2
        delete Simulation_Opt.mat
    end
    
    Day = 11;
    Exp_T = 7;
    Exp_S = 2;
    i = 2; % i is in {0,1,...,Exp_T-2}
    
    T = 2^Exp_T;
    Z21 = 2^Exp_T/4;
    Delta = (T-Z21)/(3*2^i);
    for n = 0:length(x)-1
        aux_x(1+n*Delta:(n+1)*Delta) = x(n+1);
    end
    aux_x(end+1) = aux_x(end);
    Vec = [zeros(1,Z21),aux_x,zeros(1,Z21)];
    
    fx = -1*Admissible_Solution_6('Opt',1,1,0,...
    0,0,0,Exp_S,Exp_T,1,1,...
    Vec,0,Day,i);

    gx = -1*Admissible_Solution_6('Opt',6,1,0,...
    0,0,0,Exp_S,Exp_T,1,1,...
    Vec,0,Day,i);
    
end