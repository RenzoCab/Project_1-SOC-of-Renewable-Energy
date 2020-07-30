close all
clear all
clc

Q = [1 2 3 4;
    5 6 7 8;
    9 10 11 12;
    13 14 15 16];
b = [10 20 30 40]';
d = [0.1 0.2 0.3 0.4]';
c = -10;
i = 3;
x = 2;

[newQ,newb,newc,newd] = Mat_Red(Q,b,c,d,i,x);