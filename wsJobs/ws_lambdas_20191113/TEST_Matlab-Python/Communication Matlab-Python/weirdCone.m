% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 14/11/2019

function [fx,gx] = weirdCone(x)

    disp('EVALUATING FUNCTION...');
    d  = norm(x);
    c  = floor(d+1);
    cc = floor(d);
    a = 0;
    for i = 0:cc
        a = a + i^2;
    end
    gx = 2*c*x;
    fx = c*norm(x)^2 - a;
    
end