close all
clear all
clc

for m=3:7
    for n=1:m
        
        prob(m,n) = comb(m,n)*(comb(100-m,20-n))/comb(100,20)*100;
        
    end
end
prob