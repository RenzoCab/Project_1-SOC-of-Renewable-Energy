function [X,forNeuralNetwoek] = Quad_FMC_Lambda_8(Q,b,c,d,dd,options,x0,limsA,InfTermicas)

    f = @(x) x'*dd*x + x'*d;
    
    function [cc,ceq] = nonlcon(x)
        ceq = [0];
        cc = [-x'*Q*x - x'*b + c];
%         ceq = [-x'*Q*x - x'*b + c];
%         cc = [0];
    end

    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    
    weAlterateControls = 0;
    % This flag allows us to limit the maximum power of some generator.
    if weAlterateControls
        ub(1) = 1; % Bonete T.
        ub(2) = 1; % Bonete S.
        ub(3) = 0; % Baygorria T.
        ub(4) = 0; % Baygorria S.
        ub(5) = 0; % Palmar T.
        ub(6) = 1; % SG T.
        ub(7) = 0; % Battery.
        ub(8) = 1; % Mototes Batlle.
        ub(9) = 1; % PTA.
        ub(10) = 1; % PTB.
        ub(11) = 1; % CTR.
        ub(12) = 1; % Failure.
    end
    
    if ~isempty(limsA)
        lb(7) = limsA(1);
        ub(7) = limsA(2);
    else
        lb(7) = 0;
        ub(7) = 0;
    end
    lb(8:12) = InfTermicas;

    [X,~,exitflag,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    errorCounter = 0;
    
    while exitflag ~= 1
        if exitflag == 0
            if rand(1) > 0.5
                options.Algorithm = 'interior-point';
            else
                options.Algorithm = 'sqp';
            end
        end
        if errorCounter == 1000
            error('Maximum iterations reached.');
        end
        x0 = rand(length(Q),1);
        x0(7) = 0; % The battery always in 0.
        [X,~,exitflag,~] = fmincon(f,x0,[],[],[],[],lb,ub,@nonlcon,options);
    end
        
    % NEW!!!:
    if exitflag == 1
        vectorQ = [Q(1,1),Q(1,2),Q(2,1),Q(3,3),Q(3,4),Q(4,3),Q(5,5),Q(6,6)];
        forNeuralNetwoek = {b,c,d,vectorQ,dot(X,d)};
    else
        forNeuralNetwoek = {};
    end
    
end