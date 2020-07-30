clear all;
close all;
clc;

% Day = '20180302';
% save('DayLambda.mat','Day');

load('setBadDays.mat');

for i = 1:150

    close all;
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
    
    if ~any(strcmp(setBadDays,Day))
        finiteDifferences(auxDay);
    end

end