function [X] = Opt_T_m(Q,b,c,d,T_lambda,m)

% WE DONR USE C

% Ex shared:
%     Q = [2 1/2 0;
%         1/2 1 1/4;
%         0 1/4 1];
%     b = [3 4 -1]';
%     d = [4 3 2]';
%     c = -5;
%     T_lambda = 1/5;
    
    if m == 0
        x = T_Solver(Q,b,d,T_lambda);
        if not(sum(x>1)+sum(x<0)) && isreal(x)
            x1 = x;
            return;
        else
            X{1} = 0;
            disp('No admissible solution found!');
            return;
        end
    end
    
    v=1:length(Q);
    
    X = {};
    X{1} = 0;
    [sets,bins] = Set_Bin(v,m);
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