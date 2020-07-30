close all
clear all
clc

load('WP01012017.mat');
theta0 = 0.2;
theta1 = 0.03;
alpha = 0.14;
theta = @(t) theta0*exp(-theta1*t);
t0 = 0;
T = 23;
N = 23;
dt = (T-t0)/N;
t = 0:dt:23;
lines = 3000;
Y = zeros(lines,length(t));
GP = 1383; % MW.
P = WP01012017/GP;
Y(:,1) = P(1);
xlim([t0,T]);
grid on;
hold on;

for j=1:lines
    for i=2:length(t)

        dW = sqrt(dt)*randn;
        Y(j,i) = Y(i-1) - theta(t(i-1))*(Y(i-1)-P(i-1))*dt + sqrt(2*theta(t(i-1))*alpha*Y(i-1)*(1-Y(i-1)))*dW;

    end
    plot(t,Y(j,:));
end

for i=1:length(t)
    
    E(i) = mean(Y(:,i));
    
end

p1 = plot(t,E,'r--o');
p1.LineWidth = 3;
p2 = plot(t,P,'bl--o');
p2.LineWidth = 3;

for i=1:length(t)
    
    V(i) = sum((E(i)-Y(:,i)).^2)/(lines-1);
    
end

sigma = sqrt(V);

p3 = plot(t,E+sigma,'g--o');
p3.LineWidth = 3;
p4 = plot(t,E-sigma,'g--o');
p4.LineWidth = 3;

Wind_E_sigma = [E;sigma];
save('Wind_E_sigma.mat','Wind_E_sigma');