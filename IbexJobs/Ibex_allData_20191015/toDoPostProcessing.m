close all;
clear all;
clc;

Day = '20180302';
save('Day.mat','Day');
pause(1);
errors = {};

parfor i = 1:400
    try
        pause(rand*10);
        iterationsAndOptimization;
    catch err
        errors{i} = err;
    end
end