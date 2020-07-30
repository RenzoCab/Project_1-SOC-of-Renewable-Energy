function [x1] = Quad_FMC(Q,b,c,d)

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
        c = g(x');
        ceq = c;
    end

    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    x0 = rand(1,length(Q));

    options = optimoptions('fmincon','Display','notify');
    [x1,~,flagFmincon,outputFmincon] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    fprintf('g(x1)=%f\n',g(x1'));
    aux = sprintf('%d ', x1);
    fprintf('x1=[%s]\n',aux);
        
end
