clear all;
close all;
clc;

set(0,'DefaultFigureVisible','on')
% set(0,'DefaultFigureVisible','off')

for i = 17:17

%     if isfile(['Simulation_Day_',num2str(i),'.mat'])
%         A(i) = Admissible_Solution_6(['Day_NB_',num2str(i)],4,1,0,1,1,0,2,8,1,1,0,0,i);
%     end
% 
%     if isfile(['Simulation_Day_NB_',num2str(i),'.mat'])
%         B(i) = Admissible_Solution_6(['Day_NB_',num2str(i)],2,1,0,0,0,0,2,8,1,1,0,0,i);
%         C(i) = Admissible_Solution_6(['Day_NB_',num2str(i)],4,1,0,1,1,0,2,8,1,1,0,0,i,1000);
%     end
    
%     if isfile(['Simulation_Day_',num2str(i),'.mat'])
%         D(i) = Admissible_Solution_6(['Day_',num2str(i)],2,1,1,1,1,1,2,11,1,1,0,0,i,1000);
        D(i) = Admissible_Solution_6(['Day_',num2str(i)],2,1,1,1,1,1,2,11,1,1,0,0,i,1000);
%     end
% 
%     if isfile(['Simulation_Day_NB_',num2str(i),'.mat'])
%         E(i) = Admissible_Solution_6(['Day_NB_',num2str(i)],5,1,0,0,0,0,2,8,1,1,0,0,i);
%     end
%     D = D / 100000;
%     E = E/100000;

end

return;

%% To check the lambdas:

% for j = 1:10
for j = 2:2
    for i = -200:10:200

            delete(['Test_',num2str(i),'.mat']);
            A(j,i/10+21) = Admissible_Solution_7(['Test_',num2str(i)],1,1,0,0,0,0,2,7,1,1,ones(5,1)*i/1000000,0,j,1);
            B(j,i/10+21) = Admissible_Solution_7(['Test_',num2str(i)],2,1,0,0,0,0,2,7,1,1,ones(5,1)*i/1000000,0,j,1);

    end
end

axes = [-200:10:200]/1000000;
if 0 == exist('Duality_Gaps','dir')
	mkdir('Duality_Gaps');
end

for i=2:2
figure;
plot(axes,A(i,:));
hold on;
grid minor;
plot(axes,B(i,:));
title(['Duality gap Day ',num2str(i)]);
xlim([min(axes) max(axes)]);
ylabel('Hundred Thousand Dollars');
xlabel('Constante Lambda Value');
saveas(gcf,[pwd '/Duality_Gaps/Day_',num2str(i)],'epsc');
legend('Dual Function','Primal Function');
end

for i=2:2
figure;
plot(axes,A(i,:));
hold on;
grid minor;
plot(axes,ones(41,1)*B(i,21));
title(['Duality gap Day ',num2str(i)]);
xlim([min(axes) max(axes)]);
ylabel('Hundred Thousand Dollars');
xlabel('Constante Lambda Value');
saveas(gcf,[pwd '/Duality_Gaps/Day_',num2str(i)],'epsc');
legend('Dual Function','Primal Function');
end

for i=2:2
figure;
hold on;
grid minor;
plot(axes,log10(abs((ones(1,41)*B(i,21)-A(i,:)))./(ones(1,41)*B(i,21))));
title(['Relative Difference ',num2str(i)]);
xlim([min(axes) max(axes)]);
ylabel('Relative Difference');
xlabel('Constante Lambda Value');
saveas(gcf,[pwd '/Duality_Gaps/Day_',num2str(i)],'epsc');
legend('Relative Difference');
end

for i=2:2
figure;
hold on;
grid minor;
plot(axes,(abs((ones(1,41)*B(i,21)-A(i,:)))./(ones(1,41)*B(i,21))));
title(['Duality Gap Day ',num2str(i)]);
xlim([min(axes) max(axes)]);
ylabel('Relative Difference');
xlabel('Constante Lambda Value');
saveas(gcf,[pwd '/Duality_Gaps/Day_',num2str(i)],'epsc');
legend('Guality Gap');
end
