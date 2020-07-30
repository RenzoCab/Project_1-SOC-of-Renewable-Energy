close all;
clear all;
clc;

expansion = 3;
NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

Lambda = @(A) A*ones(1,length(time));
A_TF = 0.5; % Average of previous tirbined flow.
Fcost = @(A) Z*A*A_TF;
i = 0;

% (WeUseLambda,WithWind,SaveFig,SaveU,ShowFig,ComputePlots,Use_Fmincon,expansion,Fcost,Lag,Hat_Lag).
Cost_Ini(0,0,0,0,1,1,0,expansion,Fcost(i),Lambda(i),[Lambda(i),zeros(1,Z41)]);