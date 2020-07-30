close all;
clear all;
clc;

Day = '20181201';
save('Day.mat','Day');
pause(1);

errors = {};

parfor i = 1:110
    try
        pause(rand*10);
        iterationsAndOptimization;
    end
end