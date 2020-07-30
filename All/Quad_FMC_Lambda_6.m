function [X,lambda] = Quad_FMC_Lambda_6(Q,b,c,d,dd,k,options,PhiV2Fixed,PhiV2,x0,limsA,InfTermicas)

    f = @(x) x'*dd*x + x'*d;
    
    if k(1) < 0 && abs(k(1)) > PhiV2
        k(1) = 0;
    end
    
    % Ineq:
    function [cc,ceq] = nonlcon(x)
        if PhiV2Fixed == 0
            ceq = [0,0,0];
            cc = [x(4)*k(3) + x(3)*k(2) - k(1) - 1, -(x(4)*k(3) + x(3)*k(2) - k(1)), -x'*Q*x - x'*b + c];
        elseif PhiV2Fixed == 1
            ceq = x(4)*k(3) + x(3)*k(2) - k(1) - PhiV2;
            cc = -x'*Q*x - x'*b + c;
        else
            disp('Choose PhiV2Fixed 0 or 1.');
            return;
        end
    end

    % Eq:
%     function [cc,ceq] = nonlcon(x)
%         if PhiV2Fixed == 0
%             ceq = [-x'*Q*x - x'*b + c,0];
%             cc = [x(4)*k(3) + x(3)*k(2) - k(1) - 1, -(x(4)*k(3) + x(3)*k(2) - k(1))];
%         elseif PhiV2Fixed == 1
%             ceq = [-x'*Q*x - x'*b + c, x(4)*k(3) + x(3)*k(2) - k(1) - PhiV2];
%             cc = [0,0];
%         else
%             disp('Choose PhiV2Fixed 0 or 1.');
%             return;
%         end
%     end
    
    
    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    
    if ~isempty(limsA)
        lb(7) = limsA(1);
        ub(7) = limsA(2);
    else
        lb(7) = 0;
        ub(7) = 0;
    end
    lb(8:12) = InfTermicas;

    [X,~,exitflag,~,lambda] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    
    while exitflag ~= 1
        if exitflag == 0
            options.Algorithm = 'interior-point';
        end
        x0 = rand(length(Q),1);
        x0(7) = 0; % The battery always in 0.
        [X,~,exitflag,~,lambda] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    end
    
    X(end+1) = k(3)*X(4) + k(2)*X(3) - k(1);
    
end