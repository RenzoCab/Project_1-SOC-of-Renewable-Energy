close all
clear all
clc

Tot = 29;
Num = 4:2:Tot*2+2;
for i=1:Tot
    load(['u_',num2str(3+2*i),'x',num2str(3+2*i),'.mat']);
    U{i} = u;
end

load('u_99x99.mat');
U{end+1} = u;
load('u_101x101.mat');
U{end+1} = u;
Num(end+1) = 98;
Num(end+1) = 100;

for i=1:length(U)
    for j=1:7
    N = length(squeeze(U{i}(1,:,1,4)));
    y(i,j) = U{i}(1,(N-1)/2+1,(N-1)/2+1,j);
    end
end

loglog(Num(1:end-1),abs(y(1:end-1,4)-y(end,4))/abs(y(end,4)),'-*');
title('u at the initial point');
ylabel('US$');
grid on;
xlabel('Number of partitions in Space');

%% Time:

close all
clear all
clc

load('u_21x21_T12.mat');
U{1} = u;
load('u_21x21_T24.mat');
U{2} = u;
load('u_21x21_T48.mat');
U{3} = u;
load('u_21x21_T96.mat');
U{4} = u;

for i=1:length(U)
    for j=1:7
    N = length(squeeze(U{i}(1,:,1,4)));
    y(i,j) = U{i}(1,(N-1)/2+1,(N-1)/2+1,j);
    end
end

Num = [12,24,48,96];

loglog(Num(1:end-1),abs(y(1:end-1,4)-y(end,4))/abs(y(end,4)),'-*');
title('u at the initial point');
ylabel('US$');
grid on;
xlabel('Number of partitions in Time');