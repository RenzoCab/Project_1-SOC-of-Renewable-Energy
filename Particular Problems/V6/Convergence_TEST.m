close all
clear all
clc

%% Run Test:

for i = 1:5
    
    U(i) = Admissible_Solution_4(0,0,0,2+i,3+i,0,0,0,0);
    
end

% save('U_Disconnected.mat','U'); Uncomment only to save again.
