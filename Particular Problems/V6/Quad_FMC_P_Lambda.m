function [X] = Quad_FMC_P_Lambda(Q,b,c,d,k,options)

    function y = g(x)
        y = x'*Q*x + x'*b + c;
    end
    f = @(x) d'*x';
    function [c,ceq] = nonlcon(x)
        
        ceq(1) = g(x');
        c(1) = 0;

        % ==================== Baygorria Pass Dam with VC ====================>>
        c(2) = x(5)*k(3) + x(4)*k(2) - k(1) - 1;
        ceq(2) = 0;
        c(3) = -(x(5)*k(3) + x(4)*k(2) - k(1));
        ceq(3) = 0; 
        
        % ==================== Baygorria Pass Dam without VC ====================>>
        
%         ceq(2) = x(5)*k(3) + x(4)*k(2) - k(1);
%         c(2) = 0;
        
    end

    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    x0 = ones(1,length(Q))/2;

%     options = optimoptions('fmincon','Display','notify');
    [X,~,~,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    X(8) = k(3)*X(5) + k(2)*X(4) - k(1);
    X = X';
        
end
