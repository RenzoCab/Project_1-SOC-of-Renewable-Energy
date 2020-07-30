function [] = NeededPlots(m,n,cz)

    if cz == 1
        Screensize = get(0,'screensize');
        vert = Screensize(4)/4;
        horiz = Screensize(3)/5;
        set(groot,'defaultfigureposition',[0 0 horiz vert])
    end

    for i = m:n
        figure(i);
    end

end