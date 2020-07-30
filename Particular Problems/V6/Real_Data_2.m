function [] = Real_Data_2(SaveFig)

    All = xlsread([pwd '/Data/All_Dams_JenFeb.xlsx'],1);

    for i = 1:18
        figure('Position', [10 10 900 420]);
        time = 0:1/1416:1;
        plot(time,All(:,i));
        xlim([0 1]);
        grid minor;
        xlabel('Time');
        switch i
            case 1
                ylabel('m');
                title('Up-Level Bonete');
            case 2
                ylabel('m');
                title('Down-Level Bonete');
            case 3
                ylabel('m');
                title('Delta-Level Bonete');
            case 4
                ylabel('m');
                title('Up-Level Baygorria');
            case 5
                ylabel('m');
                title('Down-Level Baygorria');
            case 6
                ylabel('m');
                title('Delta-Level Baygorria');
            case 7
                ylabel('m');
                title('Up-Level Palmar');
            case 8
                ylabel('m');
                title('Down-Level Palmar');
            case 9
                ylabel('m');
                title('Delta-Level Palmar');
            case 10
                ylabel('m');
                title('Up-Level Salto Grande');
            case 11
                ylabel('m');
                title('Down-Level Salto Grande');
            case 12
                ylabel('m');
                title('Delta-Level Salto Grande');
            case 13
                ylabel('MW');
                title('CTR');
            case 14
                ylabel('MW');
                title('PTA');
            case 15
                ylabel('MW');
                title('PTB');
            case 16
                ylabel('MW');
                title('Motores Batlle');
            case 17
                ylabel('m');
                title('Delta Bonete-Down and Baygorria-Up');
                hold on;
                plot(time,mean(All(:,i))*ones(1,length(time)),'r');
            case 18
                ylabel('m');
                title('Delta Baygrria-Down and Palmar-Up');
                hold on;
                plot(time,mean(All(:,i))*ones(1,length(time)),'r');
        end
        if 0 == exist('Dams_Plots','dir')
            mkdir('Dams_Plots');
        end
        if SaveFig == 1
            saveas(gcf,[pwd '/Dams_Plots/',num2str(i)],'epsc');
        end
    end
end

