function [X] = Quad_FMC_Lambda_8(Q,b,c,d,dd,options,x0,limsA,InfTermicas)

    f = @(x) x'*dd*x + x'*d;
    
    function [cc,ceq] = nonlcon(x)
        ceq = [0];
        cc = [-x'*Q*x - x'*b + c];
    end

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

    [X,~,exitflag,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    
    while exitflag ~= 1
        if exitflag == 0
            if rand(1) > 0.5
                options.Algorithm = 'interior-point';
            else
                options.Algorithm = 'sqp';
            end
        end
        x0 = rand(length(Q),1);
        x0(7) = 0; % The battery always in 0.
        [X,~,exitflag,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    end
        
end