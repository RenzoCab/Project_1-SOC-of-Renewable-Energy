close all;
clear all;
clc;

% Created 17/10/2019 - RC.

% DATA 2019:

bonete = readtable('Valores_de_Agua_ADME (2019).xlsx','Sheet','Gabriel Terra (Bonete)');
baygorria = readtable('Valores_de_Agua_ADME (2019).xlsx','Sheet','Rincon de Baygorria');
palmar = readtable('Valores_de_Agua_ADME (2019).xlsx','Sheet','Constitucion (Palmar)');
sg = readtable('Valores_de_Agua_ADME (2019).xlsx','Sheet','Salto Grande Uruguay');

bonete = table2array(bonete(:,2));
baygorria = table2array(baygorria(:,2));
palmar = table2array(palmar(:,2));
sg = table2array(sg(:,2));

if 1 == 0
    figure;
    histogram(bonete,100);
    title('Bonete water value histogram 2019');
    grid minor;
    xlabel('$USD/MWh$','interpreter','latex');
    saveas(gcf,[pwd '/waterValue_2019/histBonete'],'epsc');
    figure;
    histogram(baygorria,100);
    title('Baygorria water value histogram 2019');
    grid minor;
    xlabel('$USD/MWh$','interpreter','latex');
    saveas(gcf,[pwd '/waterValue_2019/histBaygorria'],'epsc');
    figure;
    histogram(palmar,100);
    title('Palmar water value histogram 2019');
    grid minor;
    xlabel('$USD/MWh$','interpreter','latex');
    saveas(gcf,[pwd '/waterValue_2019/histPalmar'],'epsc');
    figure;
    histogram(sg,100);
    title('SG water value histogram 2019');
    grid minor;
    xlabel('$USD/MWh$','interpreter','latex');
    saveas(gcf,[pwd '/waterValue_2019/histSG'],'epsc');
end

boneteNoZeros = bonete;
boneteNoZeros(boneteNoZeros == 0) = [];
baygorriaNoZeros = baygorria;
baygorriaNoZeros(baygorriaNoZeros == 0) = [];
palmarNoZeros = palmar;
palmarNoZeros(palmarNoZeros == 0) = [];
sgNoZeros = sg;
sgNoZeros(sgNoZeros == 0) = [];

disp(['Bonete nonzero-mean, median and nonzero-mode: ',num2str(mean(boneteNoZeros)),'   ',num2str(median(bonete)),'   ',num2str(mode(boneteNoZeros))]);
disp(['Baygorria nonzero-mean, median and nonzero-mode: ',num2str(mean(baygorriaNoZeros)),'   ',num2str(median(baygorria)),'   ',num2str(mode(baygorriaNoZeros))]);
disp(['Palmar nonzero-mean, median and nonzero-mode: ',num2str(mean(palmarNoZeros)),'   ',num2str(median(palmar)),'   ',num2str(mode(palmarNoZeros))]);
disp(['SG nonzero-mean, median and nonzero-mode: ',num2str(mean(sgNoZeros)),'   ',num2str(median(sg)),'   ',num2str(mode(sgNoZeros))]);

% ALL DATA:

bonete = readtable('Valores_de_Agua_ADME (PLOTS).xlsx','Sheet','Gabriel Terra (Bonete)');
baygorria = readtable('Valores_de_Agua_ADME (PLOTS).xlsx','Sheet','Rincon de Baygorria');
palmar = readtable('Valores_de_Agua_ADME (PLOTS).xlsx','Sheet','Constitucion (Palmar)');
sg = readtable('Valores_de_Agua_ADME (PLOTS).xlsx','Sheet','Salto Grande Uruguay');

bonete = table2array(bonete(:,2));
baygorria = table2array(baygorria(:,2));
palmar = table2array(palmar(:,2));
sg = table2array(sg(:,2));

boneteNoZeros = bonete;
boneteNoZeros(boneteNoZeros == 0) = [];
baygorriaNoZeros = baygorria;
baygorriaNoZeros(baygorriaNoZeros == 0) = [];
palmarNoZeros = palmar;
palmarNoZeros(palmarNoZeros == 0) = [];
sgNoZeros = sg;
sgNoZeros(sgNoZeros == 0) = [];

disp(['Bonete mean, median and nonzero-mode: ',num2str(mean(bonete)),'   ',num2str(median(bonete)),'   ',num2str(mode(boneteNoZeros))]);
disp(['Baygorria mean, median and nonzero-mode: ',num2str(mean(baygorria)),'   ',num2str(median(baygorria)),'   ',num2str(mode(baygorriaNoZeros))]);
disp(['Palmar mean, median and nonzero-mode: ',num2str(mean(palmar)),'   ',num2str(median(palmar)),'   ',num2str(mode(palmarNoZeros))]);
disp(['SG mean, median and nonzero-mode: ',num2str(mean(sg)),'   ',num2str(median(sg)),'   ',num2str(mode(sgNoZeros))]);