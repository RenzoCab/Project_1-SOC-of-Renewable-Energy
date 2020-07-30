clear all;
close all;
clc;

pc = parcluster('local');
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

load('Day.mat');
% Now we have the variable Day which is a string 'yyyymmdd'.
auxDay = Day; % We copy it (because we will now modify Day).
matFormat = datetime(Day,'InputFormat','yyyyMMdd');
% We pass from string to MATLAB format.
matFormat = matFormat + days(1);
% We add a day;
Day = datestr(matFormat,'yyyymmdd');
% We create a variable Day but with the next day.
save('Day.mat','Day');

optimumLambdaGradient(auxDay);
