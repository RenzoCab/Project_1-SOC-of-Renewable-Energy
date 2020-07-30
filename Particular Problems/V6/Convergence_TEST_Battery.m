close all
clear all
clc

%% Run Test:

for i = 2:5
    
    U(i,1) = Admissible_Solution_5(0,0,0,0,i,4+i,1,1,1,1);
    U(i,2) = Admissible_Solution_5(1,0,0,0,i,4+i,1,1,1,1);
    
end

% save('U_Battery.mat','U'); Uncomment only to save again.

%% Load Test:

load('U_Battery.mat');

%% Plot Test:

N = [2^2,2^3,2^4,2^5];

figure;
title('Normalized Cost Function');
xlabel('i'); ylabel('Normalized Cost');
hold on; grid on;
plot([2:5],U(:,1)/1e5,'*--');
plot([2:5],U(:,2)/1e5,'*--');
legend('Initial Cost without Battery','Initial Cost with Battery','location','southeast');

Y = abs(U(1:end-1,1)-U(end,1))/U(end,1);
Z = abs(U(1:end-1,2)-U(end,2))/U(end,2);

figure;
loglog(N(1:end-1),Y,'*--');
title('Convergence Test without Battery');
xlabel('Space Discretization'); ylabel('Relative Error');
grid on;
legend('Initial Cost without Battery');

figure;
loglog(N(1:end-1),Z,'*--');
title('Convergence Test with Battery');
xlabel('Space Discretization'); ylabel('Relative Error');
grid on;
legend('Initial Cost with Battery');