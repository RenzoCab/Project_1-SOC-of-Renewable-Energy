function [X] = Quad_FMC_P(Q,b,c,d,k)

% Ex shared:
%     Q = [0 0 0 0;
%         0 2 1/2 0;
%         0 1/2 1 1/4;
%         0 0 1/4 1];
%     b = [1 3 4 -1]';
%     d = [5 4 3 2]';
%     c = -5;

    function y = g(x)
        y = x'*Q*x + x'*b + c;
    end
    f = @(x) d'*x';
    function [c,ceq] = nonlcon(x)
        c(1) = g(x');
        ceq(1) = c(1);
        c(2) = x(5)*k(3) + x(4)*k(2) - k(1) - 1;
        ceq(2) = 0;
        c(3) = -(x(5)*k(3) + x(4)*k(2) - k(1));
        ceq(3) = 0;
%         c(2) = x(10)+k(1,1)-x(4)*k(1,2)-x(5)*k(1,3);
%         ceq(2) = c(2);
%         c(3) = x(11)+k(2,1)-x(6)*k(2,2)-x(7)*k(2,3);
%         ceq(3) = c(3);
    end

    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    x0 = rand(1,length(Q));

    options = optimoptions('fmincon','Display','notify');
    [X,~,flagFmincon,outputFmincon] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    fprintf('g(x1)=%f\n',nonlcon(X));
    aux = sprintf('%d ', X);
    fprintf('Min_Fmincon=[%s]\n',aux);
        
end
