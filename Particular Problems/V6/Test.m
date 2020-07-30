close all;
clear all;
clc;

Level1 = [16:0.25:28];
Power1 = [90000 90000 90000 90000 98000 98000 98000 98000 108000 108000 108000 108000 114000 114000 114000 114000 ...
    120000 120000 120000 120000 126000 126000 126000 126000 132000 132000 132000 132000 138000 138000 138000 138000 ...
    142000 142000 142000 142000 146000 146000 146000 146000 150000 150000 150000 150000 152712 152712 152712 152712 152712];
Flow1 = [654 644 634 624 692 662 652 642 703 692 681 670 699 688 679 669 696 687 ...
    677 668 706 697 688 679 692 684 676 668 691 683 676 666 677 670 664 656 667 661 655 649 660 654 648 642 648 642 636 620 624];
Level2 = [8.5:0.5:17.5];
Power2 = [60000 60000 60000 62000 64000 75000 80000 85000 89000 94000 98000 105000 106000 106000 105000 ...
    102000 100000 98000 96000];
Flow2 = [744 732 713 711 688 836 831 825 834 857 870 901 866 820 786 723 687 656 630];
Level3 = [16:31];
Power3 = [180000 200000 225000 240000 255000 285000 310000 325000 330000 330000 330000 330000 ...
    330000 330000 330000 330000];
Flow3 = [1362 1422 1519 1520 1522 1640 1707 1703 1632 1542 1466 1400 1344 1293 1247 1206];


H1 = @(V1) (-3.74)*(V1)^2 + (16.7)*(V1) + 67.7; 
H2 = 54;
H3 = @(V3) (-3.76)*(V3)^2 + (16.9)*(V3) + 26.8;
H4 = @(V4) (-19.8)*(V4)^2 + (51.5)*(V4) + 3.79;

h04 = 5.1;
h03 = 10;

FlowMax1 = @(V1) (H1(V1)-H2 <= 21)*(12.4*(H1(V1)-H2) + 447) + (H1(V1)-H2 > 21)*(-11.2*(H1(V1)-H2) + 941);
FlowMax2 = @(V3) (H2-H3(V3) <= 14)*(37.5*(H2-H3(V3)) + 376) + (H2-H3(V3) > 14)*(-81.3*(H2-H3(V3)) + 2039);
FlowMax3 = @(V3) (H3(V3)-h03 <= 21)*(46.5*(H3(V3)-h03) + 638) + (H3(V3)-h03 > 21)*(-68.1*(H3(V3)-h03) + 3272);

FlowMax4 = 4200;
MaxSpill4 = 12600;

cT1 = 1.02;
cS1 = 3.9;
d1 = @(tur,spill) cT1*tur+cS1*spill;

cT2 = 0.57;
cS2 = 1.8;
d2 = @(tur,spill) cT2*tur+cS2*spill;

cT3 = 3.77;
cS3 = 7.65;
d3 = @(tur,spill) cT3*tur+cS3*spill;

cT4 = 5.34;
cS4 = cT4*MaxSpill4/FlowMax4;
d4 = @(tur,spill) cT4*tur+cS4*spill;

eta = 9;

Power_H1 = @(x1,x2,V1) eta*x1*FlowMax1(V1)*(H1(V1)-d1(x1,x2)-H2);
Power_H2 = @(x3,x4,V3) eta*x3*FlowMax2(V3)*(H2-d2(x3,x4)-H3(V3));
Power_H3 = @(x5,V3) eta*x5*FlowMax3(V3)*(H3(V3)-d3(x5,0)-h03);
Power_H4 = @(x6,V4) eta*x6*FlowMax4*(H4(V4)-d4(x6,0)-h04);

V = [0:0.01:1];
for i = 1:length(V)
    PH1(i) = Power_H1(1,0,V(i));
    H11(i) = H1(V(i));
    PH2(i) = Power_H2(1,0,V(i));
    PH3(i) = Power_H3(1,V(i));
    H33(i) = H3(V(i));
    PH4(i) = Power_H4(1,V(i));
    H44(i) = H4(V(i));
end

if 0 == exist('Power_Plots','dir')
    mkdir('Power_Plots');
end

figure;
plot(V,PH1); grid minor;
hold on; xlim([0.2 1]);
line([0.2,1],[152712,152712]);
title('Max. Power Bonete Vs. Volume');
saveas(gcf,[pwd '/Power_Plots/1'],'epsc');

figure;
plot(V,PH2); grid minor;
hold on; xlim([0.5 1]);
line([0.5,1],[106000,106000]);
title('Max. Power Baygorria Vs. Volume');
saveas(gcf,[pwd '/Power_Plots/2'],'epsc');

figure;
plot(V,PH3); grid minor;
hold on; xlim([0.5 1]);
line([0.5,1],[330000,330000]);
title('Max. Power Palmar Vs. Volume');
saveas(gcf,[pwd '/Power_Plots/3'],'epsc');

figure;
plot(V,PH4); grid minor;
hold on; xlim([0.6 1]);
line([0.6,1],[1890000,1890000]/2);
title('Max. Power Salto Grande Vs. Volume');
saveas(gcf,[pwd '/Power_Plots/4'],'epsc');

figure;
plot(H11,PH1); grid minor;
hold on; xlim([min(H11),max(H11)]);
line([min(H11),max(H11)],[152712,152712]);
title('Max. Power Bonete Vs. Volume');
plot(54+Level1,Power1,'*');
saveas(gcf,[pwd '/Power_Plots/5'],'epsc');

figure;
plot(H33,PH2); grid minor;
hold on; xlim([min(H33),max(H33)]);
line([min(H33),max(H33)],[106000,106000]);
title('Max. Power Baygorria Vs. Volume');
plot(22.5+Level2,Power2,'*');
saveas(gcf,[pwd '/Power_Plots/6'],'epsc');

figure;
plot(H33,PH3*1.1); grid minor;
hold on; xlim([min(H33),max(H33)]);
line([min(H33),max(H33)],[330000,330000]);
title('Max. Power Palmar Vs. Volume');
plot(h03+Level3,Power3,'*');
saveas(gcf,[pwd '/Power_Plots/7'],'epsc');

figure;
plot(H44,PH4); grid minor;
hold on; xlim([min(H44),max(H44)]);
line([min(H44),max(H44)],[1890000,1890000]/2);
title('Max. Power Salto Grande Vs. Volume');
saveas(gcf,[pwd '/Power_Plots/8'],'epsc');

figure;
plot(Level1,Flow1/max(Flow1));
hold on;
plot(Level1,Power1/max(Power1));
figure;
plot(Flow1,Power1);

figure;
plot(Level2,Flow2/max(Flow2));
hold on;
plot(Level2,Power2/max(Power2));
figure;
plot(Flow2,Power2);

figure;
plot(Level3,Flow3/max(Flow3));
hold on;
plot(Level3,Power3/max(Power3));
figure;
plot(Flow3,Power3);