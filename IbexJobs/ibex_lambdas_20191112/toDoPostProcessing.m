function [] = toDoPostProcessing(WhatToDo)

    % 1 = HJB.
    % 2 = PostProcessing.

    % Day = '20180302';
    % save('Day.mat','Day');
    % save('Day_PostProcess.mat','Day');
    % save('Day_RealControls.mat','Day');

    if WhatToDo == 1
        for i=1:250
            iterationsAndOptimization(1);
            disp(i);
        end
    end

    if WhatToDo == 2
        parfor i=1:600
            try
                pause(rand*15);
                iterationsAndOptimization(2);
                disp(i);
            end
        end
    end
    
    if WhatToDo == 3
        parfor i=1:600
            try
                pause(rand*15);
                iterationsAndOptimization(3);
                disp(i);
            end
        end
    end

end