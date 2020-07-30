function [Min] = Min(Q,b,c,d,k)

    Sens_Lambda = 0.001; %sensibility of the paramether lambda
    
    Q2 = Q(2:end,2:end);
    b2 = b(2:end);
    c2 = c;
    d2 = d(2:end);
    
    C{1} = @(x) - k(1) + k(2)*x(4) + k(3)*x(5);
    
    Done = 0;
    m = 0;
    v = [1:8];
    
    while Done == 0
        
        if m == 0
            T_Lambda = b(1)/d(1);
            X = T_Solver(Q2,b2,d2,T_Lambda); % me da el vector de controles sin el fuel
            if not(sum(X>1)+sum(X<0)) && isreal(X)
                Min = Check_T(Q1,b1,c1,d1,X,C); %chekea lo virtual y lo agrega a los controles al final, tambien agrega al principio el fuel
                if isreal(Min)
                    Done = 1;
                end
            end
        end

%         if m == 1
%             
%             [sets,bins] = Set_Bin(v,m);
%             X = App_Sets(sets,bins);
%             
%         end
            
        if m > 3
            X = {};
            X{1} = 0;
            [sets,bins] = Set_Bin(v,m);
            for i = 1:length(sets(:,1))
                for j = 1:length(bins(:,1))
                    if sets(i,end) == 8
%                         [newQ,newb,newd] = App_Cond(Q1,b1,d1,i,k,fr);
                        [nQ,nb,nc,nd] = App_Cond(Q,b,c,d,bins(i,end),k);
                        sets_aux = sets(i,1:end-1);
                        bins_aux = bins(i,1:end-1);
                    else
                        nQ = Q;
                        nb = b;
                        nc = c;
                        nd = d;
                        sets_aux = sets(i,:);
                        bins_aux = bins(i,:);
                    end
                    if sets_aux(1) == 1 % this means phi_F=0, which implies NT solutions
                        nQ = nQ(2:end);
                        nb = nb(2:end);
                        nc = nc;
                        nd = nd(2:end);
                        sets_aux = sets_aux(2:end) - 1;
                        bins_aux = bins_aux(2:end);
                        if ismember(1,sets_aux) && not(ismember(2,sets_aux))
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            if ismember(3,sets_aux) && not(ismember(4,sets_aux))
                                if abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                else
                                    NoSolution = 1;
                                end
                            else
                                T_Lambda = nb(1)/nd(1);
                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,1,0); % this removes phi_S1 del sistema
                                Y = T_Solver_Bins(nQ,nb,nd,T_Lambda,sets_aux,bins_aux);
                                if isreal(Y)
                                    Y = [0,Y]; % we add phi_F=0
                                    [SQ,Sb,Sc,Sd] = Mat_Red(Q1,b1,c1,d1,3,0);
                                    Phi_S1 = (Y'*SQ*Y+Y'*Sb+Sc)/(-Y(2)*SQ(2,3)*2);
%                                     Phi_S2 = (Y'*SQ*Y+Y'*Sb+Sc)/(-Y(4)*SQ(4,5)*2);
                                    Y = [Y(1:2),Phi_S1,Y(3:end)];
                                    if C{1}(Y) < 1 && C{1}(Y) > 0 && Phi_S1 < 1 && Phi_S1 > 0
                                        X{end+1} = Y;
                                    end
                                end
                            end
                        elseif ismember(3,sets_aux) && not(ismember(4,sets_aux))
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            if ismember(1,sets_aux) && ismember(2,sets_aux)
                                T_Lambda = nb(1)/nd(1);
                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,1,0); % this removes phi_S2 del sistema
                                Y = T_Solver_Bins(nQ,nb,nd,T_Lambda,sets_aux,bins_aux);
                            elseif not(ismember(1,sets_aux)) && not(ismember(2,sets_aux))
                                T_Lambda = nb(3)/nd(3);
                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,3,0); % this removes phi_S2 del sistema
                                Y = T_Solver_Bins(nQ,nb,nd,T_Lambda,sets_aux,bins_aux);
                            else % if not(ismember(1,sets_aux)) && ismember(2,sets_aux)
                                T_Lambda = nb(2)/nd(2);
                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,2,0); % this removes phi_S2 del sistema
                                Y = T_Solver_Bins(nQ,nb,nd,T_Lambda,sets_aux,bins_aux);
                            end
                            if isreal(Y)
                                Y = [0,Y]; % we add phi_F=0
                                [SQ,Sb,Sc,Sd] = Mat_Red(Q1,b1,c1,d1,3,0);
                                Phi_S2 = (Y'*SQ*Y+Y'*Sb+Sc)/(-Y(4)*SQ(4,5)*2);
                                Y = [Y(1:4),Phi_S1,Y(5:end)];
                                if C{1}(Y) < 1 && C{1}(Y) > 0 && Phi_S2 < 1 && Phi_S2 > 0
                                    X{end+1} = Y;
                                end
                            end
                        else
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            Y = NT_Solver_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            Y = [0,Y]; % we add phi_F=0
                            if C{1}(Y) < 1 && C{1}(Y) > 0
                                X{end+1} = Y;
                            end
                        end
                    else % here Phi_F not 0
                        if  ismember(2,sets_aux) && not(ismember(3,sets_aux))
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            if  ismember(4,sets_aux) && not(ismember(5,sets_aux))
                                if abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) + abs((nb(1)/nd(1))/(nb(3)/nd(3))-1) < 2*Sens_Lambda
                                    UseMatlabMin = 1;
                                else
                                    NoSolution = 1;
                                end
                            else
                                if abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                else
                                    NoSolution = 1;
                                end
                            end
                        elseif ismember(4,sets_aux) && not(ismember(5,sets_aux))
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            if ismember(2,sets_aux) && ismember(3,sets_aux)
                                if abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                else
                                    NoSolution = 1;
                                end
                            elseif not(ismember(2,sets_aux)) && not(ismember(3,sets_aux))
                                if abs((nb(1)/nd(1))/(nb(4)/nd(4))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                else
                                    NoSolution = 1;
                                end
                            else % if not(ismember(2,sets_aux)) && ismember(3,sets_aux)
                                if abs((nb(1)/nd(1))/(nb(3)/nd(3))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                else
                                    NoSolution = 1;
                                end
                            end
                        else
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            T_Lambda = nb(1)/nd(1);
                            Y = T_Solver_Bins(nQ,nb,nd,T_Lambda,sets_aux,bins_aux);
                            if C{1}(Y) < 1 && C{1}(Y) > 0
                                X{end+1} = Y;
                            end
                        end
                    end
                end
            end
        end
        
        if length(X) > 1
            Done = 1;
        end
        
        m = m + 1;

    end
end