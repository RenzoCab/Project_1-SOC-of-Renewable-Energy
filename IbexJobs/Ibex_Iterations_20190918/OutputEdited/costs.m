function [successful] = costs()

    load('costs.mat');
    
    successful = 0;
    
    for i = 1:length(values)
                
        if length(values{i}) == 3
            
            successful = successful + 1;
            real(successful) = values{i}(1) / 1e6;
            optimalLambda(successful) = values{i}(2) / 1e6;
            lambdaZero(successful) = values{i}(3) / 1e6;
            
            if successful == 1
                realAccum(successful) = real(successful);
                optimalLambdaAccum(successful) = optimalLambda(successful);
                lambdaZeroAccum(successful) = lambdaZero(successful);
            else
                realAccum(successful) = realAccum(successful-1) + real(successful);
                optimalLambdaAccum(successful) = optimalLambdaAccum(successful-1) + optimalLambda(successful);
                lambdaZeroAccum(successful) = lambdaZeroAccum(successful-1) + lambdaZero(successful);
            end
            
            if optimalLambda(successful) > real(successful)
                successful = successful - 1;
                disp(i); % I can see which days not to add in the latex.
            end
            
        end
                    
    end
    
    figure; hold on;
    days = linspace(0,1,successful);
    P = plot(days,real,'o-'); P.LineWidth = 1;
    P = plot(days,optimalLambda,'o-'); P.LineWidth = 1;
    ylabel('Hundred Thousand USD');
    xlabel('Computed Days (Normalized)');
    legend('Historical Cost','Optimal Cost');
    title('Daily Cost Comparison');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    grid minor;
    
    figure; hold on;
    P = plot(days,realAccum,'o-'); P.LineWidth = 1;
    P = plot(days,optimalLambdaAccum,'o-'); P.LineWidth = 1;
    ylabel('Hundred Thousand USD');
    xlabel('Computed Days (Normalized)');
    legend('Historical Accumulated Cost','Optimal Accumulated Cost');
    title('Accumulated Cost Comparison');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    grid minor;

end