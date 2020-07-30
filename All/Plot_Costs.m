close all
clear all
clc

Data = xlsread('Costs 12_04_2019.xls',1);

Data(:,15) = 0;
Data(:,21) = 0;
Data(:,22) = 0;
Data(:,25) = 0;
Data(:,27) = 0;
Data(:,29) = 0;
Data(:,33) = 0;
Data(:,35) = 0;
Data(:,36) = 0;
Data(:,39) = 0;
Data(:,49) = 0;
Data(:,59) = 0;

for i = 0:7
    
    figure;
    b = bar(Data(:,1+i*7:7+i*7)');
    if i == 2 || i == 3
        legend('HJB Optimal Cost with Battery','Optimal Path with Battery',...
            'HJB Optimal Cost without Battery','Optimal Path without Battery',...
            'Historical Path','location','northwest');
    else
        legend('HJB Optimal Cost with Battery','Optimal Path with Battery',...
        'HJB Optimal Cost without Battery','Optimal Path without Battery',...
        'Historical Path');
    end
    grid minor;
    xlabel('Days'); ylabel('Hundred Thousand Dollars');
    title(['Acumulated cost during Week ',num2str(i+1)]);
%     b(1).FaceColor = [0 100/255 0];
%     b(2).FaceColor = [144/255 238/255 144/255];
%     b(3).FaceColor = [0 0 139/255];
%     b(4).FaceColor = [173/255 216/255 230/255];
%     b(5).FaceColor = [1 0 0];
    saveas(gcf,[pwd '/Historical/Cost_Week_',num2str(i+1)],'epsc');
    
end

figure;
bar(Data(:,57:59)');
legend('HJB Optimal Cost with Battery','Optimal Path with Battery',...
        'HJB Optimal Cost without Battery','Optimal Path without Battery',...
        'Historical Path');
grid minor;
xlabel('Days'); ylabel('Hundred Thousand Dollars');
title('Week 9');
saveas(gcf,[pwd '/Historical/Cost_Week_9'],'epsc');

%% Section 2:

Sav = zeros(1,59);
for i = 1:59
    if Data(4,i) ~= 0
        Sav(i) = (Data(2,i)<Data(4,i))*abs(Data(2,i)-Data(4,i))./Data(4,i);
    end
end

figure;
plot(Sav);
ylim([0,0.35]);
xlim([1,59]); grid minor;
title('Relative difference using the Battery (OP)');
xlabel('Days');
ylabel('Relative difference');
saveas(gcf,[pwd '/Historical/Relative_BatNoBat'],'epsc');

Sav2 = zeros(1,59);
for i = 1:59
    if Data(3,i) ~= 0
        Sav2(i) = (Data(1,i)<Data(3,i))*abs(Data(1,i)-Data(3,i))./Data(3,i);
    end
end

figure;
plot(Sav2);
ylim([0,0.45]);
xlim([1,59]); grid minor;
title('Relative difference using the Battery (HJB)');
xlabel('Days');
ylabel('Relative difference');
saveas(gcf,[pwd '/Historical/Relative_BatNoBat_HJB'],'epsc');

B = [sum(Data(1,:)),sum(Data(2,:)),sum(Data(3,:)),sum(Data(5,:)),sum(Data(4,:))];