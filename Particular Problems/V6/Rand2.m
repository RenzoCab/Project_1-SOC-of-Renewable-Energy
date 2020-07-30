function [R] = Rand2(a,b)

    if a == b
        R = a;
    elseif a > b
        c = a-b;
        R = b+c*rand(1);
    else
        c = b-a;
        R = a+c*rand(1);
    end

end