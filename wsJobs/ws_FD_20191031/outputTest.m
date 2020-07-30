close all;
clear all;
clc;

r{1} = testOut(1); % Only saves 'a'.

function [a,b] = testOut(z)

    a = z;
    b(1) = z;
    b(2) = z + 1;

end