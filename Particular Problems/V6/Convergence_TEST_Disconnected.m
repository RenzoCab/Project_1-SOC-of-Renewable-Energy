close all
clear all
clc

%% Run Test:

for i = 1:5
    
    U(i) = Admissible_Solution_4(0,0,0,2+i,3+i,0,0,0,0);
    
end

% save('U_Disconnected.mat','U'); Uncomment only to save again.

%% Load Test:

load('U_Disconnected.mat');

%% Plot Test:

N = [2^3,2^4,2^5,2^6,2^7];

figure;
title('Normalized Cost Function');
xlabel('i'); ylabel('Normalized Cost');
hold on; grid on;
plot([1:5],U/1e5,'*--');

Y = abs(U(1:end-1)-U(end))/U(end);

figure;
loglog(N(1:end-1),Y,'*--');
title('Convergence Test');
xlabel('Space Discretization'); ylabel('Relative Error');
grid on;