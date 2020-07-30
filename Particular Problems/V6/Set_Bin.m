function [sets,bins] = Set_Bin(v,m)

    sets = combnk(v,m);
    bins = dec2bin(2^m-1:-1:0)-'0';
    
end