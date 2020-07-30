function [X] = Quad_Solver(iQ,ib,ic,id,ik)

    function [y,grady] = quadobj(x,Q,f,c)
        y = x'*Q*x + f'*x + c;
        if nargout > 1
            grady = 2*Q*x + f;
        end
    end

    H = {-iQ,iQ,zeros(length(iQ)),zeros(length(iQ))};
    k = {-ib,ib,[0 0 0 -ik(2) -ik(3) 0 0]',[0 0 0 ik(2) ik(3) 0 0]'};
    d = {-ic,ic,ik(1),-1-ik(1)};

    function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
        jj = length(H); % jj is the number of inequality constraints
        y = zeros(1,jj);
        for i = 1:jj
            y(i) = x'*H{i}*x + k{i}'*x + d{i};
        end
        yeq = [];

        if nargout > 2
            grady = zeros(length(x),jj);
            for i = 1:jj
                grady(:,i) = 2*H{i}*x + k{i};
            end
        end
        gradyeq = [];
    end

    function hess = quadhess(x,lambda,Q,H)
        hess = Q;
        jj = length(H); % jj is the number of inequality constraints
        for i = 1:jj
            hess = hess + lambda.ineqnonlin(i)*H{i};
        end
    end

    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(x,lambda,0,H),'Display','notify');
    
    lb = zeros(1,length(iQ));
    ub = ones(1,length(iQ));

    fun = @(x)quadobj(x,0,id,0);
    nonlconstr = @(x)quadconstr(x,H,k,d);
    x0 = zeros(length(iQ),1); % Column vector.
%     [X,fval,eflag,output,lambda] = fmincon(fun,x0,...
%         [],[],[],[],lb,ub,nonlconstr,options);
    [X,~,~,~,~] = fmincon(fun,x0,...
        [],[],[],[],lb,ub,nonlconstr,options);
        
end
