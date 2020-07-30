close all;
clear all;
clc;

set(0,'DefaultFigureVisible','off');
% set(0,'DefaultFigureVisible','on');

casos = 1;
% 1 is for HJB.
% 2 is for cost comparison.

% Day = '20181201'; % We run this 2 lines to reset the day.
% save('Day.mat','Day');

if casos == 1
%     pc = parcluster('local');
%     parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end

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

if casos == 1
    try
        Optimal_Solution(auxDay,1,1,0,...
        1,1,1,3,8,1,1,...
        0,0,auxDay,999);
        Optimal_Solution(auxDay,6,1,0,...
        1,1,1,3,8,1,1,...
        0,0,auxDay,999);
    end
elseif casos == 0
    Optimal_Solution(auxDay,6,1,0,...
    1,1,1,3,8,1,1,...
    0,0,auxDay,999);
end