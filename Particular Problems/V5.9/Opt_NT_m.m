function [X] = Opt_NT_m(Q,b,c,d,m)

% Ex shared:
%     Q = [2 1/2 0;
%         1/2 1 1/4;
%         0 1/4 1];
%     b = [3 4 -1]';
%     d = [4 3 2]';
%     c = -5;

% Ex:
%     Q = [1 0 0 1/4;
%         0 2 1/2 0;
%         0 1/2 1 0;
%         1/4 0 0 1];
%     b = [1 3 4 -1]';
%     d = [5 4 3 2]';
%     c = -5;

    X = {};
    X{1} = 0;
    if m == 0
        x = NT_Solver(Q,b,c,d);
        if not(sum(x>1)+sum(x<0)) && isreal(x)
            X{end+1} = x;
            return;
        else
            disp('No admissible solution found!');
            return;
        end
    end
    
    v=1:length(Q);
    function y = g(x)
        y = x'*Q*x + x'*b + c;
    end
    f = @(x) d'*x';
    function [c,ceq] = nonlcon(x)
        c = g(x');
        ceq = c;
    end
    
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
            y = NT_Solver(newQ,newb',newc,newd');               
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