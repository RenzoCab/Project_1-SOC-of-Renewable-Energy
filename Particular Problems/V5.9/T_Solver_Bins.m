function [X] = T_Solver_Bins(nQ,nb,nd,T_Lambda,sets,bins)

    Y = T_Solver(nQ,nb,nd,T_Lambda);
    if isreal(Y) == 1
        cont1 = 1;
        cont2 = 1;
        for l=1:length(Q)+length(bins)
            if cont1<=length(bins) && l==sets(cont1)
                X(l) = bins(cont1);
                cont1 = cont1 +1;
            else
                X(l) = Y(cont2);
                cont2 = cont2 +1;
            end
        end
        if not(sum(X>1)+sum(X<0))
            X = 1i;
        end
    else
        X = 1i;
    end
    
end