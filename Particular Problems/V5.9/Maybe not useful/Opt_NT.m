function [x1] = Opt_NT(Q,b,c,d,UseFmincon)

    m = 1;
    Done = 0;

%     Q = [1 1/2 0;
%         1/2 1 0;
%         0 0 1];
    % b = [1 0 1]';
    % d = [1 1 1]';
    % c = -1;

    v=1:length(Q);

    function y = g(x)
        y = x'*Q*x + x'*b + c;
    end
    f = @(x) d'*x;
    function [c,ceq] = nonlcon(x)
        c = g(x');
        ceq = c;
    end

    x = NT_Solver(Q,b,c,d);
    if not(sum(x>1)+sum(x<0))% && abs(g(x))<1e-9
        Done = 1;
    end

    while Done == 0 || m > length(v)
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
                y = NT_Solver(newQ,newb',newc,newd');               
                if isreal(y) == 1
                    cont1 = 1;
                    cont2 = 1;
                    for l=1:length(Q)
                        if cont1<=m && l==sets(i,cont1)
                        %if ismember(l,sets(i,:))
                            x(l) = bins(j,cont1);
                            cont1 = cont1 +1;
                        else
                            x(l) = y(cont2);
                            cont2 = cont2 +1;
                        end
                    end
                    if not(sum(x>1)+sum(x<0)) %&& abs(g(x))<1e-9
                        Done = 1;
                        X{end+1} = x;
                    end
                end    
            end
        end
        m = m + 1;
    end

    if Done == 1
        CostVal = zeros(length(X)-1,1);
        for i=2:length(X)
            CostVal(i-1) = f(X{i});
        end
        ind_op = find(CostVal == max(CostVal)); % Here we may use min.
        x1 = X{ind_op+1};
    elseif m > length(v)
        disp('No admissible solution found!');
        return;
    end

    fprintf('g(x1)=%f\n',g(x1));
    aux = sprintf('%d ', x1);
    fprintf('x1=[%s]\n',aux);

    % ==================== Verification (fmincon) ====================>>>

    if UseFmincon == 1
        h = @(x) sum(-x);
        lb = zeros(1,length(Q));
        ub = ones(1,length(Q));
        x0 = rand(1,length(Q));
        
        options = optimoptions('fmincon','Display','notify');
        [x2,~,flagFmincon,outputFmincon] = fmincon(h,x0,[],[],[],[],lb,ub,@nonlcon,options);
        fprintf('g(x2)=%f\n',g(x2'));
        aux = sprintf('%d ', x2);
        fprintf('x2=[%s]\n',aux);
    elseif UseFmincon ~= 0
        disp('Choose UseFmincon between 0 or 1.');
        return;
    end
    x2'

end