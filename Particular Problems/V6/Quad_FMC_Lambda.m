function [X] = Quad_FMC_Lambda(Q,b,c,d,k,options,PhiV2Fixed,PhiV2,x0)

    f = @(x) x'*d;
    
    function [cc,ceq] = nonlcon(x)
        if PhiV2Fixed == 0
            ceq = [x'*Q*x + x'*b + c,0];
            cc = [x(5)*k(3) + x(4)*k(2) - k(1) - 1, -(x(5)*k(3) + x(4)*k(2) - k(1))];
        elseif PhiV2Fixed == 1
            ceq = [x'*Q*x + x'*b + c, x(5)*k(3) + x(4)*k(2) - k(1) - PhiV2];
            cc = [0,0];
        else
            disp('Choose PhiV2Fixed 0 or 1.');
            return;
        end
    end

    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    
    [X,~,exitflag,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    
    while exitflag ~= 1
        x0 = zeros(7,1);
        [X,~,exitflag,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    end
    
    X(8) = k(3)*X(5) + k(2)*X(4) - k(1);
    
end
