function [Min] = Minimization(Q,b,d,c,k)

    Q1 = Q(2:end,2:end);
    b1 = b(2:end);
    d1 = d(2:end);
    T_Lambda = b(1)/d(1);
    Done = 0;
    m = 0;
    C{1} = @(x) - k(1,1) + k(1,2)*x(4) + k(1,3)*x(5);
    C{2} = @(x) - k(2,1) + k(2,2)*x(6) + k(2,3)*x(7);
    k10 = [k(1,1)/k(1,3) (1+k(1,1))/k(1,3); k(2,1)/k(2,3) (1+k(2,1))/k(2,3)];
    fr = [k(1,2)/k(1,3) k(2,2)/k(2,3)];
    
    while Done == 0

        if m == 0
            X = Opt_T_m(Q1,b1,c,d1,T_Lambda,m);
            if length(X) > 1
                Min = Check_T(Q,b,c,d,X(2:end),C);
                if isreal(Min)
                    Done = 1;
                end
            end
        end
        
        if m == 1
            Y{1} = Opt_NT_m(Q1,b1,c,d1,m-1);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            Y{2} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            Y{3} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,5,k10(2,1),fr(2));
            Y{4} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,5,k10(2,2),fr(2));
            Y{5} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            
            Y{6} = Opt_T_m(Q1,b1,c,d1,T_Lambda,m);
            
            X = {};
            X{1} = 0;
            
            for i = 1:length(Y)
                if length(Y{i}) > 1
                    X = {X,Y{2:end}};
                end
            end
            
            if length(X) > 1
                Min = Check_T(Q,b,c,d,X(2:end),C);
                if isreal(Min)
                    Done = 1;
                end
            end
        end
        
        if m == 2
            Y = {};
            Y{1} = Opt_NT_m(Q1,b1,c,d1,m-1);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            Y{2} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{3} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            Y{4} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{5} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,5,k10(2,1),fr(2));
            Y{6} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{7} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,5,k10(2,2),fr(2));
            Y{8} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{9} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,1),fr(2));
            Y{10} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,2),fr(2));
            Y{11} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,1),fr(2));
            Y{12} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,2),fr(2));
            Y{13} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);

            Y{14} = Opt_T_m(Q1,b1,c,d1,T_Lambda,m);
            
            X = {};
            X{1} = 0;
            
            for i = 1:length(Y)
                if length(Y{i}) > 1
                    X = {X,Y{2:end}};
                end
            end
            
            if length(X) > 1
                Min = Check_T(Q,b,c,d,X(2:end),C);
                if isreal(Min)
                    Done = 1;
                end
            end
        end
        
        if m > 2
            Y = {};
            Y{1} = Opt_NT_m(Q1,b1,c,d1,m-1);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            Y{2} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{3} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            Y{4} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{5} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,5,k10(2,1),fr(2));
            Y{6} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{7} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,5,k10(2,2),fr(2));
            Y{8} = Opt_T_m(Q2,b2,c,d2,T_Lambda,m-1);
            Y{9} = Opt_NT_m(Q2,b2,c,d2,m-2);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,1),fr(2));
            Y{10} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            Y{11} = Opt_NT_m(Q3,b3,c,d3,m-3);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,1),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,2),fr(2));
            Y{12} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            Y{13} = Opt_NT_m(Q3,b3,c,d3,m-3);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,1),fr(2));
            Y{14} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            Y{15} = Opt_NT_m(Q3,b3,c,d3,m-3);
            
            [Q2,b2,d2] = App_Cond(Q1,b1,d1,3,k10(1,2),fr(1));
            [Q3,b3,d3] = App_Cond(Q2,b2,d2,4,k10(2,2),fr(2));
            Y{16} = Opt_T_m(Q3,b3,c,d3,T_Lambda,m-2);
            Y{17} = Opt_NT_m(Q3,b3,c,d3,m-3);

            Y{18} = Opt_T_m(Q1,b1,c,d1,T_Lambda,m);
            
            X = {};
            X{1} = 0;
            
            for i = 1:length(Y)
                if length(Y{i}) > 1
                    X = {X,Y{2:end}};
                end
            end
            
            if length(X) > 1
                Min = Check_T(Q,b,c,d,X(2:end),C);
                if isreal(Min)
                    Done = 1;
                end
            end
        end
        
        m = m + 1;
        
    end

end