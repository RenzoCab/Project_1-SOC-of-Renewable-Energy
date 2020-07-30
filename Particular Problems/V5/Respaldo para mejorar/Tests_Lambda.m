%% TEST 0: A*ones(1,time):

close all
clear all
clc

expansion = 2;
Int = 0.005;
A = linspace(-Int,Int,201);

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

Lambda = @(A) A*[zeros(1,Z41),ones(1,length(time)-Z41)];

WithOut_Lambda = Cost_Ini(0,0,0,0,0,0,0,expansion,Lambda(A(1)),[Lambda(A(1)),zeros(1,Z41)]);

parfor i=1:length(A)
    U(i) = Cost_Ini(1,0,0,0,0,0,0,expansion,Lambda(A(i)),[Lambda(A(i)),zeros(1,Z41)]);
    Y(i) = WithOut_Lambda;
end

close all;
set(groot,'defaultFigureVisible','on');
plot(A,U);
hold on;
plot(A,Y);
grid on;
title('Dual cost function');
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
MinD = -1e-3;
MaxD = 1e-3;
A = linspace(MinD,MaxD,25);
B = linspace(MinD,MaxD,25);

NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

Lambda = @(A,B) (A + B*sin(3*pi*time)).*[zeros(1,Z41),ones(1,length(time)-Z41)];

WithOut_Lambda = Cost_Ini(0,0,0,0,0,0,0,expansion,Lambda(A(1),B(1)),[Lambda(A(1),B(1)),zeros(1,Z41)]);

for i=1:length(A)
    parfor j=1:length(B)
    
        U(i,j) = Cost_Ini(1,0,0,0,0,0,0,expansion,Lambda(A(i),B(j)),[Lambda(A(i),B(j)),zeros(1,Z41)]);
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
save('Exp_1','U1');

saveas(gcf,'Cost_Lambda_1','epsc');

%% TEST 2: A0 + Ai*sin(2*pi*i*time) + Bi*cos(2*pi*i*time) i=1,2,...,10:

close all
clear all
clc

global expansionSin

expansionSin = 2;

NT = 4*2^expansionSin; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 21; % Number of coefficient.

WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansionSin,zeros(NT+1,1),[zeros(1,NT+1),zeros(1,Z41)]);

options = optimset('PlotFcns',@optimplotfval,'TolFun',1.0000e-7);
set(groot,'defaultFigureVisible','on');
figure(1000);
Lambda_Op = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansionSin,LambdaFSinCos(x),[LambdaFSinCos(x),zeros(1,Z41)]),[0.1,zeros(n,1)'],options);
grid on;
saveas(gcf,'Cost_Lambda_2_OP','epsc');

set(0,'CurrentFigure',1000);
grid on;
title('Optimal Lambda function');
xlabel('Time (t)');
saveas(gcf,'Cost_Lambda_2','epsc');
U2 = [WithOut_Lambda,Lambda_Op];

save('Exp_2','U2');

%% TEST 3: Polynomilas:

close all
clear all
clc

global expansionPol

expansionPol = 2;

NT = 4*2^expansionPol; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 21; % Number of coefficient.

WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansionPol,zeros(NT+1,1),[zeros(1,NT+1),zeros(1,Z41)]);

options = optimset('PlotFcns',@optimplotfval,'TolFun',1.0000e-7);
set(groot,'defaultFigureVisible','on');
figure(1000);
Lambda_Op = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansionPol,LambdaFPol(x),[LambdaFPol(x),zeros(1,Z41)]),ones(n,1)*0,options);
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

global expansionF
expansionF = 2;

NT = 4*2^expansionF; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 4*2^expansionF+1; % Number of coefficient.
WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansionF,zeros(NT+1,1),[zeros(1,NT+1),zeros(1,Z41)]);

set(groot,'defaultFigureVisible','on');
figure(1000);

options = optimset('PlotFcns',@optimplotfval,'TolFun',1.e-12,'MaxFunEvals',1e6,'MaxIter',1e6);
options.StepTolerance = 1e-16;
[Lambda_Op_1,~] = fminunc(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansionF,LambdaF(x),[LambdaF(x),zeros(1,Z41)]),ones(n,1)*0.1,options);

%====

% Int = 20;
% A = linspace(-Int,Int,1000);
% 
% Lambda = @(A) [Lambda_Op_1(1:3)',A,Lambda_Op_1(5:end)'];
% 
% parfor i=1:length(A)
%     U(i) = Cost_Ini(1,0,0,0,0,0,0,expansionF,Lambda(A(i)),[Lambda(A(i)),zeros(1,Z41)]);
%     Y(i) = WithOut_Lambda;
% end
% 
% close all;
% set(groot,'defaultFigureVisible','on');
% plot(A,U);
% hold on;
% plot(A,Y);
% grid on;
% title('Cost function at initial value');
% xlabel('A');
% ylabel('M USD');
% legend('Relaxed','No Relaxed')
% [maxU,idxU] = max(U(:));
% vline(A(idxU),'r',num2str(A(idxU)));
% xlim([-Int,Int]);

%====

[Lambda_Op_2,~] = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansionF,LambdaF(x),[LambdaF(x),zeros(1,Z41)]),(rand(n,1)-0.5)*0.1,options);

grid on;
saveas(gcf,'Cost_Lambda_4_OP','epsc');

set(0,'CurrentFigure',1000);
grid on;
title('Optimal Lambda function');
xlabel('Time (t)');
saveas(gcf,'Cost_Lambda_4','epsc');
U4 = [WithOut_Lambda,Lambda_Op_1',Lambda_Op_2'];

save('Exp_4','U4');

%% Step by step:

close all
clear all
clc

global expansionFin
global num
num = 1;
expansionFin = 1;

NT = 4*2^expansionFin; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.
n = 4*2^expansionFin+1-Z41; % Number of coefficient.
WithOut_Lambda = -Cost_Ini(0,0,0,0,0,0,0,expansionFin,zeros(NT+1,1),[zeros(1,NT+1),zeros(1,Z41)]);

set(groot,'defaultFigureVisible','on');
figure(1000);

options = optimset('PlotFcns',@optimplotfval,'TolFun',1.e-12,'MaxFunEvals',1e6,'MaxIter',1e6);
options.StepTolerance = 1e-16;
Lambda_Op_1 = 0;
 
for i=1:6
    [Lambda_Op_1,~] = fminsearch(@(x) -Cost_Ini(1,0,0,0,0,0,0,expansionFin,LambdaFin(x),[LambdaFin(x),zeros(1,Z41)]),Lambda_Op_1,options);
    num = num + 1;
    if i~=1
        A = linspace(0,1,num-1);
        B = linspace(0,1,num);
        Lambda_Op_1 = interp1(A,Lambda_Op_1,B);
    else
        Lambda_Op_1 = ones(1,2)*Lambda_Op_1;
    end
end

grid on;
saveas(gcf,'Cost_Lambda_5_OP','epsc');

set(0,'CurrentFigure',1000);
grid on;
title('Optimal Lambda function');
xlabel('Time (t)');
saveas(gcf,'Cost_Lambda_5','epsc');
U5 = [WithOut_Lambda,Lambda_Op_1];

save('Exp_5','U5');

%% Functions:

function y = LambdaFSinCos(coeffs)

    global expansionSin
    % I define this again:
    NT = 4*2^expansionSin;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    % =====
        
    y = 0;
    for i = 1:floor(length(coeffs)/2)+1
        y = y + coeffs(i)*cos(2*pi*(i-1)*time);
    end
    for i = floor(length(coeffs)/2)+2:length(coeffs)
        y = y + coeffs(i)*sin(2*pi*(i-floor(length(coeffs)/2)+1)*time);
    end
    
    y = y.*[zeros(1,Z41),ones(1,length(time)-Z41)];
    
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.0001);
end

function y = LambdaF(coeffs)

    global expansionF
    % I define this again:
    NT = 4*2^expansionF;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    % =====
        
    temp = repmat(coeffs*1e-1,[1,ceil((NT+1)/length(coeffs))])';
    y = temp(1:NT+1);
    y = y.*[zeros(1,Z41),ones(1,length(time)-Z41)];
%     y = y.*[zeros(1,Z41),ones(1,length(time)-Z41-1),0];
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.0001);
end

function y = LambdaFPol(coeffs)

    global expansionPol
    % I define this again:
    NT = 4*2^expansionPol;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    % =====
        
    y = 0;
    for i = 1:length(coeffs)
        y = y + coeffs(i)*time.^(i-1);
    end
    y = y.*[zeros(1,Z41),ones(1,length(time)-Z41)];
%     y = y.*[zeros(1,Z41),ones(1,length(time)-Z41-1),0];
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.0001);

end

function y = LambdaFin(coeffs)

    global expansionFin num
    NT = 4*2^expansionFin;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    Cant = NT - Z41 + 1;
   
    sec = floor(Cant/num);
    endsec = Cant-sec*num;
    y = zeros(1,Cant);
    
    for i=1:num
        y(1+(i-1)*sec:i*sec) = coeffs(i)';
    end
    if endsec ~=0
        y(num*sec+1:end) = coeffs(end)';
    end
    y = [zeros(1,Z41),y];
    
    set(0,'CurrentFigure',1000);
    plot(time,y);
    pause(0.00001);
end