clear all;
close all;
clc;

Day = '20180302';
save('DayLambda.mat','Day');
save('DayLambda_1.mat','Day');
save('DayLambda_2.mat','Day');
save('DayLambda_3.mat','Day');
save('DayLambda_4.mat','Day');
save('DayLambda_5.mat','Day');
jumpDays = 35;

try
	pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end

load('setBadDays.mat');

for i = 1:10

    load('DayLambda_2.mat');
    % Now we have the variable Day which is a string 'yyyymmdd'.
    auxDay = Day; % We copy it (because we will now modify Day).
    matFormat = datetime(Day,'InputFormat','yyyyMMdd');
    % We pass from string to MATLAB format.
    matFormat = matFormat + days(jumpDays);
    % We add a day;
    Day = datestr(matFormat,'yyyymmdd');
    % We create a variable Day but with the next day.
    save('DayLambda_2.mat','Day');
    
    if ~any(strcmp(setBadDays,auxDay))
%         try
            optimumLambdaGradient(auxDay,2);
%         end
    end

end