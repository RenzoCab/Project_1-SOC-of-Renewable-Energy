function [] = Optimum_Lambda_Sections(expansion_T,expansion_S,Day)

    NT = 2^expansion_T;
    Tmax = 1; Tmin = 0;
    dt = (Tmax-Tmin)/NT;
    time = Tmin:dt:Tmax;
    T = length(time);
    t21 = 1/4; Z21 = floor(length(time)*t21);
    t32 = 10/24; Z32 = floor(length(time)*t32);
    End = 0;
    Div = 1;
    Count = 1;
    figure(50);
    Th1 = 1e-5;
    
%     function [vector] = VecLambdas(coeff)
%         vector = [];
%         CantZ21 = T - Z21;
%         dZ21 = floor((T - Z21)/Div);
%         for k = 1:Div
%             vector = [vector,ones(1,dZ21)*coeff(k)];
%             CantZ21 = CantZ21 - dZ21;
%         end
%         if CantZ21 ~= 0
%             vector = [vector,ones(1,CantZ21)*coeff(Div)];
%         end
%         vector = [zeros(1,Z21),vector(1:T-Z21),zeros(1,Z21)];
%         if dZ21 == 1
%         End = 1;
%         end
%         set(0,'CurrentFigure',50);
%         plot(1:length(vector),vector);
%         xlim([1 length(vector)]);
%         title('Optimal Lambda');
%         grid on;
%         xlabel('Normalized Time');
%         pause(0.0001);
%     end

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
        if dZ21 == 1
        End = 1;
        end
%         vector = vector*0;
%         vector(Z21+1+4+4+4+4+4+4+4+4+4+4+4+4:Z21+4+4+4+4+4+4+4+4+4+4+4+4+4) = coeff;
        set(0,'CurrentFigure',50);
        plot(1:length(vector),vector);
        xlim([1 length(vector)]);
        title('Optimal Lambda');
        grid on;
        xlabel('Normalized Time');
        pause(0.0001);
    end
    
    Op_Lambda = 0;
    DG = 10;
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-5,'TolX',1e-6,'MaxFunEvals',200,'MaxIter',200);

    while End ~= 1 && DG > 0.03 && Count ~= 2
        
        [Op_Lambda,fval] = fminsearch(@(x) -Admissible_Solution_6(['Lambda_Day_',num2str(Day)],1,1,0,0,0,0,expansion_S,expansion_T,1,1,...
            subsref(VecLambdas(x),struct('type','()','subs',{{1:T+Z21}})),...
            0,Day),Op_Lambda,options);
            
        Info{Count,1} = -fval;
        Info{Count,2} = Admissible_Solution_6(['Lambda_Day_',num2str(Day)],2,1,0,0,0,0,expansion_S,expansion_T,1,1,...
            subsref(VecLambdas(Op_Lambda),struct('type','()','subs',{{1:T+Z21}})),...
            0,Day);
        Info{Count,3} = Div;
        Info{Count,4} = Op_Lambda;
        Info{Count,5} = VecLambdas(Op_Lambda);
        DG = abs((Info{Count,2}-Info{Count,1}))/Info{Count,2};
        
        save(['Op_Lambda_2',num2str(Count)],'Info');
        
        Div = Div*2;
        Count = Count + 1;
        
        Aux = Op_Lambda;
        for i = 0:length(Aux)-1
            Op_Lambda(1+2*i:2+2*i) = Aux(i+1);
        end
        
    end
end