close all
clear all
clc

profile clear
profile on

Error = 0;
Sin_Error = -1;
options = optimoptions('fmincon','Display','notify','Algorithm','sqp');

while Error == 0 && Sin_Error < 50

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
    
    Min_fmincon_1 = Quad_FMC_Erik(Q,b,c,d,k,options);
    Min_fmincon_1(8) = [];
    Cost1 = Min_fmincon_1'*d;
    Cond1_1 = Min_fmincon_1'*Q*Min_fmincon_1 + Min_fmincon_1'*b + c;
    Cond1_2 = Min_fmincon_1(5)*k(3) + Min_fmincon_1(4)*k(2) - k(1);

    Min_mine = Min_Reverse(Q,b,c,d,k)';
    Min_mine(8) = [];
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

profile viewer
% profsave(profile('info'),'Minimizator_TEST_4')