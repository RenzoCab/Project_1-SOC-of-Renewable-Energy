function [] = iterationsAndOptimization(WhatToDo)

    set(0,'DefaultFigureVisible','off');
    % set(0,'DefaultFigureVisible','on');

    % 1 is for HJB.
    % 2 is for cost comparison.

    % Day = '20181201'; % We run this 2 lines to reset the day.
    % save('Day.mat','Day');

    if WhatToDo == 1
        try
            pc = parcluster('local');
            parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
        end
    end
    
    if WhatToDo == 1
        load('Day.mat');
    elseif WhatToDo == 2
        load('Day_PostProcess.mat');
    elseif WhatToDo == 3
        load('Day_RealControls.mat');
    end
    % Now we have the variable Day which is a string 'yyyymmdd'.
    auxDay = Day; % We copy it (because we will now modify Day).
    matFormat = datetime(Day,'InputFormat','yyyyMMdd');
    % We pass from string to MATLAB format.
    matFormat = matFormat + days(1);
    % We add a day;
    Day = datestr(matFormat,'yyyymmdd');
    % We create a variable Day but with the next day.
    if WhatToDo == 1
        save('Day.mat','Day');
    elseif WhatToDo == 2
        save('Day_PostProcess.mat','Day');
    elseif WhatToDo == 3
        save('Day_RealControls.mat','Day');
    end

    if WhatToDo == 1
        try
            Optimal_Solution(auxDay,1,1,0,...
            1,1,1,3,8,1,1,...
            0,0,auxDay,999);
            Optimal_Solution(auxDay,6,1,0,...
            1,1,1,3,8,1,1,...
            0,0,auxDay,999);
        end
    elseif WhatToDo == 2
        try
            % SavePlots = 0!!!
            Optimal_Solution(auxDay,6,1,0,...
            1,1,0,3,8,1,1,...
            0,0,auxDay,999);
        end
    elseif WhatToDo == 3
        try
            Optimal_Solution(auxDay,7,1,0,...
            1,1,1,3,8,1,1,...
            0,0,auxDay,999);
        end
    end
end