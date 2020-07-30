close all
clear all
clc

load('u_V1_4_V4_4_W_4_T_4.mat');
U{1} = u;
load('u_V1_8_V4_8_W_4_T_8.mat');
U{2} = u;
load('u_V1_16_V4_16_W_4_T_16.mat');
U{3} = u;
load('u_V1_32_V4_32_W_4_T_32.mat');
U{4} = u;
load('u_V1_64_V4_64_W_4_T_64.mat');
U{5} = u;
load('u_V1_128_V4_128_W_4_T_128.mat');
U{6} = u;
load('u_V1_256_V4_256_W_4_T_256.mat');
U{7} = u;

load('u_V1_4_V4_4_W_8_T_4.mat');
Y{1} = u;
load('u_V1_8_V4_8_W_8_T_8.mat');
Y{2} = u;
load('u_V1_16_V4_16_W_8_T_16.mat');
Y{3} = u;
load('u_V1_32_V4_32_W_8_T_32.mat');
Y{4} = u;
load('u_V1_64_V4_64_W_8_T_64.mat');
Y{5} = u;
load('u_V1_128_V4_128_W_8_T_128.mat');
Y{6} = u;
load('u_V1_256_V4_256_W_8_T_256.mat');
Y{7} = u;

T4 = [4,8,16,32,64,128,256];
T8 = [4,8,16,32,64,128,256];

for i=1:length(U)
    
    P4(i) = U{i}(1,floor((2*2^i)/2)+1,floor((2*2^i)/2)+1,3);
    
end

for i=1:length(Y)
    
    P8(i) = Y{i}(1,floor((2*2^i)/2)+1,floor((2*2^i)/2)+1,5);
    
end

loglog(T4,abs(P4-P4(end))/abs(P4(end)));
hold on;
grid on;
loglog(T8,abs(P8-P8(end))/abs(P8(end)));
title('Relative error');
xlabel('Number of discretizations in time and water');
legend('4 disc. in wind','8 disc. in wind');
saveas(gcf,'Relative_error','epsc')

figure;
semilogx(T4,P4);
grid on;
hold on;
semilogx(T8,P8);
title('Cost function at inicial point');
ylabel('US$');
xlabel('Number of discretizations in time and water');
legend({'4 disc. in wind','8 disc. in wind'},'Position',[0.75 0.8 0 0]);
saveas(gcf,'Cost_function_at_inicial_point','epsc')