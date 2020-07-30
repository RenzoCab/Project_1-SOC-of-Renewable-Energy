function [X] = Quad_FMC_Erik_FixV2(Q,b,c,d,k,options,V2)

    f = @(x) d'*x';
    function [cc,ceq] = nonlcon(x)
        ceq = [x*Q*x' + x*b + c,x(5)*k(3) + x(4)*k(2) - k(1) - V2];
        cc = [0,0];
    end

    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    x0 = ones(1,length(Q))/2;

    [X,~,~,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    X(8) = k(3)*X(5) + k(2)*X(4) - k(1);
    X = X';
        
end
