%% TEST 0: A*ones(1,time):

close all
clear all
clc

expansion = 2;
Int = 0.1;
A = linspace(-Int,Int,200);

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

Lambda = @(A) A*ones(1,length(time));
A_TF = 0.5; % Average of previous tirbined flow.
Fcost = @(A) Z*A*A_TF;

WithOut_Lambda = Cost_Ini(0,0,0,0,0,0,0,expansion,Fcost(A(1)),Lambda(A(1)),[Lambda(A(1)),zeros(1,Z41)]);

parfor i=1:length(A)
    U(i) = Cost_Ini(1,0,0,0,0,0,0,expansion,Fcost(A(i)),Lambda(A(i)),[Lambda(A(i)),zeros(1,Z41)]);
    Y(i) = WithOut_Lambda;
end

close all;
set(groot,'defaultFigureVisible','on');
plot(A,U);
hold on;
plot(A,Y);
grid on;
title('Cost function at initial value');
xlabel('A');
ylabel('M USD');
legend('Relaxed','No Relaxed')
[maxU,idxU] = max(U(:));
vline(A(idxU),'r',num2str(A(idxU)));
xlim([-Int,Int]);
U0 = {A,U,Y};
save('Exp_0','U0');

saveas(gcf,'Cost_Lambda_0','epsc');

%% TEST 1: A + B*sin(3*pi*time):

close all
clear all
clc

expansion = 2;
MinD = -1;
MaxD = 1;
A = linspace(MinD,MaxD,25);
B = linspace(MinD,MaxD,25);

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

Lambda = @(A,B) A + B*sin(3*pi*time);
A_TF = 0.5;
Fcost = @(A,B) (-B/(3*pi) + A*Z + B*cos(2*pi*Z)/(3*pi))*A_TF;

WithOut_Lambda = Cost_Ini(0,0,0,0,0,0,0,expansion,Fcost(A(1),B(1)),Lambda(A(1),B(1)),[Lambda(A(1),B(1)),zeros(1,Z41)]);

for i=1:length(A)
    parfor j=1:length(B)
    
        U(i,j) = Cost_Ini(1,0,0,0,0,0,0,expansion,Fcost(A(i),B(j)),Lambda(A(i),B(j)),[Lambda(A(i),B(j)),zeros(1,Z41)]);
        Y(i,j) = WithOut_Lambda;
        
    end
end

[a,b] = meshgrid(A,B);
close all;
set(groot,'defaultFigureVisible','on');
surf(a,b,U);
hold on;
s = surf(a,b,Y,'FaceAlpha',0.5);
s.EdgeColor = 'none';
title('Cost function at initial value');
xlabel('A');
ylabel('B');
zlabel('M USD');
xlim([MinD MaxD]);
ylim([MinD MaxD]);
view(140,20);
U1 = {a,b,U,Y};
% save('Exp_1','U1');

% saveas(gcf,'Cost_Lambda_1','epsc');

%% TEST 2: A0 + Ai*sin(2*pi*i*time) + Bi*sin(2*pi*i*time) i=1,2,...,10:

close all
clear all
clc

expansion = 1;

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 21; % Number of coefficient.

WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansion,1,zeros(n,1),[zeros(1,n),zeros(1,Z41)]);

options = optimset('PlotFcns',@optimplotfval);
set(groot,'defaultFigureVisible','on');
figure(1000);
Lambda_Op = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansion,FcostFSinCos(x),LambdaFSinCos(x),[LambdaFSinCos(x),zeros(1,Z41)]),ones(n,1)*0,options);
grid on;
saveas(gcf,'Cost_Lambda_2_OP','epsc');

set(0,'CurrentFigure',1000);
grid on;
title('Optimal Lambda function');
xlabel('Time (t)');
saveas(gcf,'Cost_Lambda_2','epsc');
U2 = [WithOut_Lambda,Lambda_Op'];

save('Exp_2','U2');

%% TEST 3: Polynomilas:

close all
clear all
clc

expansion = 1;

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 21; % Number of coefficient.

WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansion,1,zeros(n,1),[zeros(1,n),zeros(1,Z41)]);

options = optimset('PlotFcns',@optimplotfval);
set(groot,'defaultFigureVisible','on');
figure(1000);
Lambda_Op = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansion,FcostFPol(x),LambdaFPol(x),[LambdaFPol(x),zeros(1,Z41)]),ones(n,1)*0,options);
grid on;
saveas(gcf,'Cost_Lambda_3_OP','epsc');

set(0,'CurrentFigure',1000);
grid on;
title('Optimal Lambda function');
xlabel('Time (t)');
saveas(gcf,'Cost_Lambda_3','epsc');
U3 = [WithOut_Lambda,Lambda_Op'];

save('Exp_3','U3');

%% TEST 4: Deltas:

close all
clear all
clc

expansion = 1;

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 4*2^(expansion)+1; % Number of coefficient.

WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansion,1,zeros(n,1),[zeros(1,n),zeros(1,Z41)]);

options = optimset('PlotFcns',@optimplotfval,'TolFun',1.0000e-07);
set(groot,'defaultFigureVisible','on');
figure(1000);
Lambda_Op = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansion,FcostF(x),LambdaF(x),[LambdaF(x),zeros(1,Z41)]),ones(n,1)*0,options);
grid on;
saveas(gcf,'Cost_Lambda_4_OP','epsc');

set(0,'CurrentFigure',1000);
grid on;
title('Optimal Lambda function');
xlabel('Time (t)');
saveas(gcf,'Cost_Lambda_4','epsc');
U4 = [WithOut_Lambda,Lambda_Op'];

save('Exp_4','U4');

%% Functions:

function y = LambdaFSinCos(coeffs)

    % I define this again:
    expansion = 2; 
    NT = 4*2^expansion;
    time = 0:1/NT:1;
    % =====
        
    y = 0;
    for i = 1:floor(length(coeffs)/2)+1
        y = y + coeffs(i)*cos(2*pi*(i-1)*time);
    end
    for i = floor(length(coeffs)/2)+2:length(coeffs)
        y = y + coeffs(i)*sin(2*pi*(i-floor(length(coeffs)/2)+1)*time);
    end
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.001);
end

function y = FcostFSinCos(coeffs)

    % I define this again:
    expansion = 2;
    NT = 4*2^expansion; %R.
    time = 0:1/NT:1; %R.
    Z = 1/3;
    Z41 = floor(length(time)*Z);
    % =====
    
    A_TF = 0.5;
    Hat_Lambda = [LambdaFSinCos(coeffs),zeros(1,Z41)]*A_TF;
    y = sum(Hat_Lambda(2:Z41)/NT);
    
end

function y = LambdaF(coeffs)

    % I define this again:
    expansion = 1; 
    NT = 4*2^expansion;
    time = 0:1/NT:1;
    % =====
        
    temp = repmat(coeffs*1e-1,[1,ceil((NT+1)/length(coeffs))])';
    y = temp(1:NT+1);
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.001);
end

function y = FcostF(coeffs)

    % I define this again:
    expansion = 1;
    NT = 4*2^expansion; %R.
    time = 0:1/NT:1; %R.
    Z = 1/3;
    Z41 = floor(length(time)*Z);
    % =====
    
    A_TF = 0.5;
    Hat_Lambda = [LambdaF(coeffs),zeros(1,Z41)]*A_TF;
    y = sum(Hat_Lambda(2:Z41)/NT);
    
end

function y = LambdaFPol(coeffs)

    % I define this again:
    expansion = 1; 
    NT = 4*2^expansion;
    time = 0:1/NT:1;
    % =====
        
    y = 0;
    for i = 1:length(coeffs)
        y = y + coeffs(i)*time.^(i-1);
    end
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.001);

end

function y = FcostFPol(coeffs)

    % I define this again:
    expansion = 1;
    NT = 4*2^expansion; %R.
    time = 0:1/NT:1; %R.
    Z = 1/3;
    Z41 = floor(length(time)*Z);
    % =====
    
    A_TF = 0.5;
    Hat_Lambda = [LambdaFPol(coeffs),zeros(1,Z41)]*A_TF;
    y = sum(Hat_Lambda(2:Z41)/NT);
    
end