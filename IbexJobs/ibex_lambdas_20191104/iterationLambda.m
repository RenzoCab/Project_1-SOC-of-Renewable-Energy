clear all;
close all;
clc;

% Day = '20180302';
% save('DayLambda.mat','Day');

try
	pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end

load('setBadDays.mat');

for i = 1:10

    load('DayLambda.mat');
    % Now we have the variable Day which is a string 'yyyymmdd'.
    auxDay = Day; % We copy it (because we will now modify Day).
    matFormat = datetime(Day,'InputFormat','yyyyMMdd');
    % We pass from string to MATLAB format.
    matFormat = matFormat + days(1);
    % We add a day;
    Day = datestr(matFormat,'yyyymmdd');
    % We create a variable Day but with the next day.
    save('DayLambda.mat','Day');

    % To do Testing:
%     auxDay = '20180302';
    
    if ~any(strcmp(setBadDays,auxDay))
        try
            optimumLambdaGradient(auxDay,mod(i,5));
        end
    end

end