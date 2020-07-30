function [] = Optimum_Lambda(expansion)

    NT = 2^expansion;
    Tmax = 1; Tmin = 0; dt = (Tmax-Tmin)/NT; time = Tmin:dt:Tmax;
    T = length(time);
    t21 = 1/4; Z21 = floor(length(time)*t21);
    t32 = 1/6; Z32 = floor(length(time)*t32);
    End = 0;
    Div = 1;
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
        
        CantZ32 = T - Z32;
        dZ32 = floor((T - Z32)/Div);
        for k = 1:Div
            vector = [vector,ones(1,dZ32)*coeff(Div+k)];
            CantZ32 = CantZ32 - dZ32;
        end
        if CantZ32 ~= 0
            
            vector = [vector,ones(1,CantZ32)*coeff(2*Div)];
        
        end
        
        vector = [zeros(1,Z21),vector(1:T-Z21),zeros(1,Z21),zeros(1,Z32),vector(T-Z21+1:end),zeros(1,Z32)];
        
        if dZ21 == 1 || dZ32 == 1
            End = 1;
        end
        
        set(0,'CurrentFigure',50);
        plot(1:length(vector),vector);
        xlim([1 length(vector)]);
        title('Optimal Lambda');
        grid on;
        xlabel('Normalized Time');
        pause(0.0001);

    end
    
    load('Op_Lambda_3.mat'); % EXTRA.
    C = [0,0];
    C = Info{3}; % EXTRA.
    Div = Info{2} + 1; % EXTRA.
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1,'TolX',1e-7,'MaxFunEvals',200,'MaxIter',200);

%     while Div ~= 4
    while End ~= 1 % EXTRA.
        
        Op_Lambda = fminsearch(@(x) -Admissible_Solution_4(0,0,0,expansion-2,expansion,1,1,...
            subsref(VecLambdas(x),struct('type','()','subs',{{1:T+Z21}})),...
            subsref(VecLambdas(x),struct('type','()','subs',{{T+Z21+1:2*T+Z21+Z32}}))),...
            C,options);
        
        Info{1} = Admissible_Solution_4(0,0,0,expansion-2,expansion,1,1,...
            subsref(VecLambdas(Op_Lambda),struct('type','()','subs',{{1:T+Z21}})),...
            subsref(VecLambdas(Op_Lambda),struct('type','()','subs',{{T+Z21+1:2*T+Z21+Z32}})));
        Info{2} = Div;
        Info{3} = Op_Lambda;
        Info{4} = VecLambdas(Op_Lambda);
        
        save(['Op_Lambda_',num2str(Div)],'Info');
        
        Div = Div + 1;
        
        for j = 1:length(Op_Lambda)
            C(2*j-1:2*j) = Op_Lambda(j);
        end
        
    end
end