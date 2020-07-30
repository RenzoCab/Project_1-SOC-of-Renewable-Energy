close all
clear all
clc

load('Test.mat');


Min_fmincon = Quad_FMC_P(Q(1:end-1,1:end-1),b(1:end-1),c,d(1:end-1),k)

Min_mine = Min(Q,b,c,d,k)