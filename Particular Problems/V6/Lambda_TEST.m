close all
clear all
clc

expansion = 3;
Disc = 2*2^expansion;

Lambda21 = ones(1,Disc+1);
Lambda32 = ones(1,Disc+1);

t21 = 1/4;
Z21 = floor(length(Lambda21)*t21);
t32 = 1/6;
Z32 = floor(length(Lambda21)*t32);

Lambda21 = [Lambda21,zeros(1,Z21)];
Lambda32 = [Lambda32,zeros(1,Z32)];
Lambda21(1:Z21) = 0;
Lambda32(1:Z32) = 0;

U = Cost_Ini_Deterministic_Lambda(0,0,expansion,1,1,Lambda21,Lambda32);
U1 = Cost_Ini_Deterministic(0,0,expansion);

A = [-2:0.25:2]*0.05;
B = [-1:0.5:1]*0.001;

L = zeros(length(A),length(B));

parfor i = 1:length(A)
    P = zeros(1,length(A));
    for j = 1:length(B)
        P(j) = Cost_Ini_Deterministic_Lambda(0,0,expansion,1,1,Lambda21*A(i),Lambda32*B(j))
    end
    L(i,:) = P;
end

[vplot1,vplot3] = meshgrid(A,B);
surf(vplot1,vplot3,L')
title('Cost as function of the Lambdas');
xlabel('\lambda_{21}');
ylabel('\lambda_{32}');

saveas(gcf,['Lambda_',num2str(expansion)],'epsc');