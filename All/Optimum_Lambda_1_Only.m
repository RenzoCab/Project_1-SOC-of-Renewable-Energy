function [] = Optimum_Lambda_1_Only(expansion_T,expansion_S,Day)

    NT = 2^expansion_T;
    Tmax = 1; Tmin = 0;
    dt = (Tmax-Tmin)/NT;
    time = Tmin:dt:Tmax;
    T = length(time);
    t21 = 1/4; Z21 = floor(length(time)*t21);
    Div = T - Z21;
    Count = 1;
    figure(50);
    
    function [vector] = VecLambdas(coeff)

        vector = [];
        
        CantZ21 = T - Z21;
        dZ21 = floor((T - Z21)/Div);
        for k = 1:Div
            vector = [vector,ones(1,dZ21)*coeff(k)];
            CantZ21 = CantZ21 - dZ21;
        end
        if CantZ21 ~= 0
            vector = [vector,ones(1,CantZ21)*coeff(Div)];
        end
        
        vector = [zeros(1,Z21),vector(1:T-Z21),zeros(1,Z21)];
        
        set(0,'CurrentFigure',50);
        plot(1:length(vector),vector);
        xlim([1 length(vector)]);
        title('Optimal Lambda');
        grid on;
        xlabel('Normalized Time');
        pause(0.0001);

    end
    
    Op_Lambda = zeros(1,T - Z21);
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-3,'TolX',1e-6,'MaxFunEvals',2000,'MaxIter',2000);
        
    [Op_Lambda,fval] = fminsearch(@(x) -Admissible_Solution_6(['Lambda_Day_',num2str(Day)],1,1,0,0,0,0,expansion_S,expansion_T,1,1,...
        subsref(VecLambdas(x),struct('type','()','subs',{{1:T+Z21}})),...
        0,Day),Op_Lambda,options);

    Info{1} = -fval;
    Info{2} = Admissible_Solution_6(['Lambda_Day_',num2str(Day)],2,1,0,0,0,0,expansion_S,expansion_T,1,1,...
        subsref(VecLambdas(Op_Lambda),struct('type','()','subs',{{1:T+Z21}})),...
        0,Day);
    Info{3} = Div;
    Info{4} = Op_Lambda;
    Info{5} = VecLambdas(Op_Lambda);

    save(['Op_Lambda_1_Only',num2str(Count)],'Info');
        
end