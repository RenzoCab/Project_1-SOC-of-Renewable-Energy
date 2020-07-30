close all
clear all
clc

% profile clear
% profile on

Error = 0;
Sin_Error = -1;

time0 = 0;
time1 = 0;
time2 = 0;

while Error == 0 && Sin_Error < 200

    P = [2000,160,108,333,900]; % All the powers.
    K(1) = P(1);
    K(2:5) = P(2:5)/30;
    d = [P(1)*Rand2(10,20)*2, P(2)*Rand2(10,20), P(2)*Rand2(10,20), P(3)*Rand2(10,20), P(3)*Rand2(10,20), P(4)*Rand2(10,20), P(4)*Rand2(10,20)]'/1000;
    b = [K(1)*Rand2(0.8,1.2), K(2)*Rand2(20,30), 0, K(3)*Rand2(20,30), 0, K(4)*Rand2(20,30), K(5)*Rand2(20,30)]';
    Q = zeros(7);
    Q(2,2) = -K(2)*Rand2(1,2);
    Q(2,3) = -K(2)*Rand2(1,2)/2;
    Q(3,2) = Q(2,3);
    Q(4,4) = -K(3)*Rand2(1,2);
    Q(4,5) = -K(3)*Rand2(1,2)/2;
    Q(5,4) = Q(4,5);
    Q(6,6) = -K(4)*Rand2(1,2);
    Q(7,7) = -K(5)*Rand2(1,2);
    k = [1 1 1];
    c = -Rand2(2500,200);
    
    tic
    Min_fmincon_2 = Quad_FMC_P(Q,b,c,d,k,optimoptions('fmincon','Display','notify','Algorithm','interior-point'));
    time0 = time0 + toc;
    
    tic
    Min_fmincon_1 = Quad_FMC_P(Q,b,c,d,k,optimoptions('fmincon','Display','notify','Algorithm','sqp'));
    time1 = time1 + toc;
    Min_fmincon_1(8) = [];
    
    Cost1 = Min_fmincon_1'*d;
    Cond1_1 = Min_fmincon_1'*Q*Min_fmincon_1 + Min_fmincon_1'*b + c;
    Cond1_2 = Min_fmincon_1(5)*k(3) + Min_fmincon_1(4)*k(2) - k(1);
   
    H = {-Q,Q,zeros(length(Q)),zeros(length(Q))};
    ik = {-b,b,[0 0 0 -k(2) -k(3) 0 0]',[0 0 0 k(2) k(3) 0 0]'};
    id = {-c,c,k(1),-1-k(1)};
    
    lb = zeros(1,length(Q));
    ub = ones(1,length(Q));
    fun = @(x)quadobj(x,0,d,0);
    nonlconstr = @(x)quadconstr(x,H,ik,id);
    
    tic
    Min_mine = fmincon(fun,zeros(1,length(Q))',...
        [],[],[],[],lb,ub,nonlconstr,optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(x,lambda,0,H),'Display','notify'));
    time2 = time2 + toc;
    
    Cost2 = Min_mine'*d;
    Cond2_1 = Min_mine'*Q*Min_mine + Min_mine'*b + c;
    Cond2_2 = Min_mine(5)*k(3) + Min_mine(4)*k(2) - k(1);
    
    if Cost1 < Cost2
        
        if abs(Min_fmincon_1 - Min_mine) > 0.01
            Error = 1;
        end
        
    end
    
    Sin_Error = Sin_Error + 1;
    disp(Sin_Error);

end

% profile viewer
% profsave(profile('info'),'Minimizator_TEST_4')

function hess = quadhess(x,lambda,Q,H)
    hess = Q;
    jj = length(H); % jj is the number of inequality constraints
    for i = 1:jj
        hess = hess + lambda.ineqnonlin(i)*H{i};
    end
end

function [y,grady] = quadobj(x,Q,f,c)
    y = x'*Q*x + f'*x + c;
    if nargout > 1
        grady = 2*Q*x + f;
    end
end

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