clear all;
close all;
clc;

for i = 1:59

B = Admissible_Solution_6(['Day_',num2str(i)],3,1,1,...
    1,1,1,2,11,1,1,...
    0,0,i);

end