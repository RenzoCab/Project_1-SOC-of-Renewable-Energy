function [Val] = Int4D(vals,x,divs)
    
    dx = 1/divs;
    m = 4;
    X = ~(dec2bin(2^m-1:-1:0)-'0').*1;
    X(1,1) = 1e-6; % To remove warning.
    [ww,xx,yy,zz] = ndgrid(0:dx:1);
    xq = [ww(:) xx(:) yy(:) zz(:)];
    vq = griddatan(X,vals,xq,'linear',{'QJ'});
    vq = reshape(vq,size(xx));
    
    dom = 0:dx:1;
    for i = 1:4
        [~,index] = min(abs(dom-x(i)));
        coef(i) = index;
    end
    
    Val = vq(coef(1),coef(2),coef(3),coef(4));

end

% We are interpolating in an unitary hypercube of dimension 4, the idea is
% to give in vals the values in each note (0000,0001,...,1111), in divs the
% amout of divisionf we want in each coordinate (they are from 0 to 1, with divs=3 we would have
% 0,0.333,0.666,1) and in x the coordinates of the interpolated point (for example
% [0.5,0.5,0.5,0.5] would give us the interpolated function in the center
% of the cube).