function [X] = App_Sets(sets,bins)


        

    for i = 1:length(sets(:,1))
        for j = 1:length(bins(:,1))
            for k = 1:length(bins(1,:))
                if k == 1
                    [newQ,newb,newc,newd] = Mat_Red(Q,b,c,d,sets(i,k)+1-k,bins(j,k));
                else
                    [newQ,newb,newc,newd] = Mat_Red(newQ,newb,newc,newd,sets(i,k)+1-k,bins(j,k));
                end     
            end
            y = T_Solver(newQ,newb',newd',T_lambda);
            if isreal(y) == 1
                cont1 = 1;
                cont2 = 1;
                for l=1:length(Q)
                    if cont1<=m && l==sets(i,cont1)
                        x(l) = bins(j,cont1);
                        cont1 = cont1 +1;
                    else
                        x(l) = y(cont2);
                        cont2 = cont2 +1;
                    end
                end
                if not(sum(x>1)+sum(x<0))
                    X{end+1} = x;
                end
            end    
        end
    end




end