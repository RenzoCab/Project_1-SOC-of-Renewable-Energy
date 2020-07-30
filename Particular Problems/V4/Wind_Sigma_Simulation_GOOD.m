close all
clear all
clc

load('WP01012017.mat');

Sigma_W(1) = 0;
dt = 1; % Is 1 hour.

GW = 1383; % MW.
% p = WP01012017/GW;
p = WP01012017/GW*2*0.8+0.1;
% theta_W0 = 0.2;
theta_W0 = 0.4;
% theta_W1 = 0.03;
theta_W1 = 0.01;
% alpha_W = 0.14;
alpha_W = 0.1;
theta_W = @(t) theta_W0*exp(-theta_W1*t);

for s=2:26
    
    t = s-1;
    Sigma_W(s) = Sigma_W(s-1) + dt*(   -(2*theta_W(t-1)+2*alpha_W*theta_W(t-1)...
        *(1-p(s-1))*p(s-1))*Sigma_W(s-1)   +   2*alpha_W*theta_W(t-1)*((1-p(s-1)*p(s-1))^2)   );
    
end

Sigma_W = sqrt(Sigma_W);

plot([0:25],Sigma_W);
figure;
plot([0:25],p);
xlim([0,25]);

figure;
hold on;
xlim([0,25]);
grid on;
title('Normalized wind power and SD')
xlabel('time (hs)');
saveas(gcf,'Wind_with_SD','epsc');

% === Derivative ===

for d = 2:25
    
    dp(1) = 0;
    der_sig_W(d) = (Sigma_W(d+1)-Sigma_W(d-1))/(2*dt);
    dp(d) = (p(d+1)-p(d-1))/(2*dt);
        
end

for i = 1:20
    
    for j = 2:25
        Y(i,1) = p(1);
        Y(i,j) = Y(i,j-1) + 1*(dp(j-1)-theta_W(j-1))*(Y(i,j-1)-p(j-1))+...
            sqrt(2*theta_W(j-1)*alpha_W*Y(i,j-1)*(1-Y(i,j-1)))*randn;
    end
    
    T = real(Y(i,:));
    
    if (min(T)>=0) && (max(T)<=1)
        Y2 = interp1([0:24],T,[0:0.1:24],'cubic');
        plot([0:0.1:24],Y2);
    end
    
end
xlim([0,24]);
grid minor;
ylabel('Normalized Wind Power');
title('Wind Power Paths');
box on;
saveas(gcf,'Wind_Paths','epsc');

figure;
p2 = interp1([0:24],p(1:25),[0:0.1:24],'cubic');
P = plot([0:0.1:24],p2);
P.LineWidth = 2;
hold on;
p2 = interp1([0:24],p(1:25)+Sigma_W(1:25)',[0:0.1:24],'cubic');
P = plot([0:0.1:24],p2);
P.LineWidth = 2;
p2 = interp1([0:24],p(1:25)-Sigma_W(1:25)',[0:0.1:24],'cubic');
P = plot([0:0.1:24],p2);
P.LineWidth = 2;
xlim([0,24]);
box on;
grid minor;
xlabel('Time (h)');
ylabel('Normalized Wind Power');
title('Expected Wind Power and Standard Deviation')
legend('Expected Wind Power','Expected + Standard Deviation','Expected - Standard Deviation')
saveas(gcf,'Wind_Space','epsc');

figure;
plot([2:25],der_sig_W(2:25));
xlim([2,25]);

figure;
plot([2:25],dp(2:25));
xlim([2,25]);

Wind_Data = [p(1:25)';dp(1:25);Sigma_W(1:25);der_sig_W(1:25)];
save('Wind_Data.mat','Wind_Data');