% close all
% clear all
% clc

load('Test.mat');

Min_fmincon_1 = Quad_FMC_P(Q,b,c,d,k);

% Min_fmincon_2 = Quad_Solver(Q,b,c,d,k)

Min_mine = Min(Q,b,c,d,k);