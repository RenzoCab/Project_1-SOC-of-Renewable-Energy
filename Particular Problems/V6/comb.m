function [x] = comb(m,n)

    x = factorial(m)/(factorial(m-n)*factorial(n));

end