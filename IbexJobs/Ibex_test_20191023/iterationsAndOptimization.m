close all;
clear all;
clc;

set(0,'DefaultFigureVisible','off');
casos = 1;
auxDay = '20181201'; % We run this 2 lines to reset the day.

if casos == 1
    try
    pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
    end
end

% Optimal_Solution has FARFOR inside:
try
    Optimal_Solution(auxDay,1,1,0,...
    1,1,1,3,8,1,1,...
    0,0,auxDay,999);
end