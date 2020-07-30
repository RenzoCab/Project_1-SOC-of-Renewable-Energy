function [x1] = Check_T(Q,b,c,d,X,C)

% Ex shared:
%     C{1} = @(x) 5 - 2*x(1)^2 - x(2)^2 - x(3)^2 - x(1)*x(2) ...
%     - x(2)*x(3)/2 - 3*x(1) - 4*x(2) + x(3);
    
    f = @(x) d'*x';
    function y = g(x)
        y = x'*Q*x + x'*b + c;
    end

    for i=1:length(X)
        x = (-c - X{i}(2:end)'*Q(2:end,2:end)*X{i}(2:end) - b(2:end)*X{i}(2:end)) / b(1); % This solver the fuel and add it to the vector of controls
        X{i} = [x,X{i}];
    end
        
    Y{1} = 0;
    for i = 1:length(X)
        for j = 1:length(C)
            y(j) = C{j}(X{i}); % Here we have the condition over the virtual control.
        end
        if  not(sum(y>1)+sum(y<0))
            Y{end+1} = [X{i},y];
        end
    end
    
    if length(Y) > 1
        
        CostVal = zeros(length(Y)-1,1);
        for i=2:length(Y)
            CostVal(i-1) = f(Y{i}(1:end-length(C)));
        end
        ind_op = find(CostVal == min(CostVal));
        x1 = Y{ind_op+1};

        fprintf('g(x1)=%f\n',g(x1'));
        aux = sprintf('%d ', x1);
        fprintf('x1=[%s]\n',aux);

    else
        
        x1 = i;
%         disp('No admissible solution found!');
%         return;
        
    end

end

% Here y=C(X) must be a cobtrol, since after we check if it is between 0 amd 1