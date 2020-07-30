close all;
clear all;
clc;

try
    A = @(t) 3;
    A(b);
catch err
    errorMessage = err.message;
    disp(errorMessage);
end