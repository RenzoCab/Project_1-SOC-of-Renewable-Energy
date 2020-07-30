function [Minimum] = Min_Reverse(Q,b,c,d,k)

    Sens_Lambda = 0.001; % Sensibility of the parameter Lambda. 
    EndWhile = 0;
    UseMatlabMin = 0;
    m = 7;
    SolExists = 0;
    T_Lambda_F = b(1)/d(1);
    WeIncreased = 0;
    Sensibility = 0.01;
    v = 1:8; % El octavo es la desigualdad en el control virtual, pero no afecta Q, b, c ni d.
    
    while EndWhile == 0
        
        m = m - 1;
        X = {};
        X{1} = 0;
            
        [sets,bins] = Set_Bin(v,m);
        for i = 1:length(sets(:,1))
            for j = 1:length(bins(:,1))

%                     if (ismember(4,sets(i,:))&&(ismember(5,sets(i,:)))) || (ismember(4,sets(i,:))&&(ismember(8,sets(i,:)))) || (ismember(5,sets(i,:))&&(ismember(8,sets(i,:))))
                % We will replace the previous line by using the next
                % piece of code.
                ToCheck = 0;
                BinsToCheck = [0,0,0];
                BinsAux = 1;
                for s = 1:length(sets(i,:))
                    if sets(i,s) == 4 || sets(i,s) == 5 || sets(i,s) == 8
                        ToCheck = ToCheck + sets(i,s);
                        BinsToCheck(BinsAux) = bins(j,s);
                        BinsAux = BinsAux + 1;
                    end
                end
                if ToCheck > 8

%                         [x4,x5,x8] = Condition(sets(i,:),bins(j,:),k);
                    % We replace by the next function.
                    [x4,x5,x8] = Condition_Faster(Sensibility,k,ToCheck,BinsToCheck);

                    if x4 ~= 2
                        [nQ,nb,nc,nd] = Mat_Red(Q,b,c,d,5,x5); % First 5.
                        [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,4,x4); % After 4.
                        sets_aux = sets(i,:);
                        bins_aux = bins(j,:);

                        if sets_aux(end) == 8
                            sets_aux = sets_aux(1:end-1);
                            bins_aux = bins_aux(1:end-1);
                        end

                        if ToCheck == 9 || ToCheck == 13 || ToCheck == 17
                            bins_aux(sets_aux==5) = [];
                            sets_aux(sets_aux==5) = [];
                        end
                        if ToCheck == 9 || ToCheck == 12 || ToCheck == 17
                            bins_aux(sets_aux==4) = [];
                            sets_aux(sets_aux==4) = [];

                        end % We transform [A,4,5,B] --> [A,B].
                        for u = 1:length(sets_aux)
                            if sets_aux(u) > 3
                                sets_aux(u) = sets_aux(u) - 2;
                            end
                        end % We transform [1,2,3,6,7] --> [1,2,3,4,5].

                        % Now we have the next system without T2 and S2.

                        if ~isempty(sets_aux) && sets_aux(1) == 1 && bins_aux(1) == 0 % This means Phi_F=0 which implies NT solutions.

                            nQ = nQ(2:end,2:end);
                            nb = nb(2:end);
                            nd = nd(2:end);
                            sets_aux = sets_aux(2:end) - 1; % We remove Phi_F and [2,3,4,5] --> [1,2,3,4].
                            bins_aux = bins_aux(2:end);
                            x1 = 0;

                            if ~isempty(sets_aux) && sets_aux(1) == 1 && (length(sets_aux) == 1 || sets_aux(2) ~= 2) % This means T1 fixed and S1 free.

                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,1,bins_aux(1)); % We reduce T1.
                                x2 = bins_aux(1); % This is T1.
                                sets_aux = sets_aux(2:end);
                                bins_aux = bins_aux(2:end);
                                sets_aux = sets_aux - 1;

                                T_Lambda_S1 = nb(1)/nd(1); % Lambda for S1.
                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,1,0); % We reduce S1.
                                sets_aux = sets_aux - 1;

                                if length(nQ) ~= length(sets_aux)

                                    [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                    Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_S1,sets_aux,bins_aux);
                                    % Y should have the values of T3 and T4
                                    % which are x6 and x7.

                                    if isreal(Y)
                                        Y = [x1,x2,x4,x5,Y];
                                        [sQ,sb,sc,~] = Mat_Red(Q,b,c,d,3,0);
                                        x3 = (Y*sQ*Y'+Y*sb+sc)/(-2*Q(2,3)*x2);
                                        Y = [Y(1:2),x3,Y(3:end),x8];
                                        if not(sum(Y>1) + sum(Y<0))
                                            X{end+1} = Y;
                                        end
                                    end

                                else

                                    x2 = bins(j,sets(i,:)==2);
                                    x6 = bins(j,sets(i,:)==6);
                                    x7 = bins(j,sets(i,:)==7);
                                    Y = [x1,x2,0,x4,x5,x6,x7];
                                    x3 = (Y*Q*Y'+Y*b+c)/(-2*Q(2,3)*x2);
                                    Y(3) = x3;
                                    Y = [Y,x8];
                                    if not(sum(Y>1) + sum(Y<0))
                                        X{end+1} = Y;
                                    end

                                end

                            else % This means NO(T1 fixed and S1 free).
                                [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                Y = NT_Solver_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                % Y should have the values of T1, S1,
                                % T3 and T4.  Y = [x2,x3,x6,x7].
                                if isreal(Y)
                                    Y = [x1,Y(1:2),x4,x5,Y(3:4),x8];
                                    if not(sum(Y>1) + sum(Y<0))
                                        X{end+1} = Y;
                                    end
                                end
                            end

                        elseif ~isempty(sets_aux) && sets_aux(1) == 1 && bins_aux(1) == 1

                            % We ignore this case.

                        else % Here Phi_F is free.

                            if ~isempty(sets_aux) && sets_aux(1) == 2 && (length(sets_aux) == 1 || sets_aux(2) ~= 3)
                                % Puedo ignorar c y nc.
                                nQ = nQ(2:end,2:end);
                                nb = nb(2:end);
                                nd = nd(2:end);
                                [~,nb,~,nd] = Mat_Red(nQ,nb,nc,nd,1,bins_aux(1)); % Reduzco T1.
                                T_Lambda_S1 = nb(1)/nd(1);
                                if (T_Lambda_F/T_Lambda_S1 > 0) && abs((T_Lambda_F/T_Lambda_S1) - 1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end

                            else 

                                [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,1,0);
                                sets_aux = sets_aux - 1;

                                if length(nQ) ~= length(sets_aux)

                                    [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                    Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_F,sets_aux,bins_aux);
                                    % Y should have the values of T1, S1,
                                    % T3 and T4. Y = [x2,x3,x6,x7].

                                    if isreal(Y)
                                        Y = [0,Y(1:2),x4,x5,Y(3:4)];
                                        x1 = (Y*Q*Y'+Y*b+c)/(-b(1));
                                        Y(1) = x1;
                                        Y = [Y,x8];
                                        if not(sum(Y>1) + sum(Y<0))
                                            X{end+1} = Y;
                                        end
                                    end

                                else

                                    x2 = bins(j,sets(i,:)==2);
                                    x3 = bins(j,sets(i,:)==3);
                                    x6 = bins(j,sets(i,:)==6);
                                    x7 = bins(j,sets(i,:)==7);
                                    Y = [0,x2,x3,x4,x5,x6,x7];
                                    x1 = (Y*Q*Y'+Y*b+c)/(-b(1));
                                    Y(1) = x1;
                                    Y = [Y,x8];
                                    if not(sum(Y>1) + sum(Y<0))
                                        X{end+1} = Y;
                                    end

                                end

                            end

                        end % We finish checking Phi_F.

                    end

                    % No solution (the condition fails).

                elseif sets(i,end) == 8 % 4 and 5 are not in sets(i,:).

                    % We know that 4 is not in sets(i,:), but we do all
                    % like it was.
                    [nQ,nb,nc,nd] = App_Cond(Q,b,c,d,bins(j,end),k,4); % 4 is the position of T2.
                    x8 = bins(j,end);
                    sets_aux = sets(i,1:end-1);
                    bins_aux = bins(j,1:end-1);

                    for u = 1:length(sets_aux)
                        if sets_aux(u) > 5
                            sets_aux(u) = sets_aux(u) - 1;
                        end
                    end % We transform [1,2,3,5,6,7] --> [1,2,3,4,5,6].

                    % Now we have the next system without S2.

                    if sets_aux(1) == 1 && bins_aux(1) == 0 % This means Phi_F=0 which implies NT solutions.

                        nQ = nQ(2:end,2:end);
                        nb = nb(2:end);
                        nd = nd(2:end);
                        sets_aux = sets_aux(2:end) - 1; % We remove Phi_F and [2,3,4,5,6] --> [1,2,3,4,5].
                        bins_aux = bins_aux(2:end);
                        x1 = 0;

                        if ~isempty(sets_aux) && sets_aux(1) == 1 && (length(sets_aux) == 1 || sets_aux(2) ~= 2) % This means T1 fixed and S1 free.

                            [nQ,nb,nc,nd] = Mat_Red(Q,b,c,d,1,bins_aux(1)); % We reduce T1.
                            x2 = bins_aux(1); % This is T1.
                            sets_aux = sets_aux(2:end);
                            bins_aux = bins_aux(2:end);

                            T_Lambda_S1 = nb(1)/nd(1); % Lambda for S1.
                            [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_S1,sets_aux,bins_aux);
                            % Y should have the values of T2, T3 and T4
                            % which are x4, x6 and x7.

                            if isreal(Y)
                                x5 = (x8 - Y(1)*k(2) + k(1)) / k(3); % Y(1) = x4.
                                Y = [x1,x2,Y(1),x5,Y(2:3)]; % Y has all except x3.
                                [sQ,sb,sc,~] = Mat_Red(Q,b,c,d,3,0);
                                x3 = (Y*sQ*Y'+Y*sb+sc)/(-x2*Q(2,3)*2);
                                Y = [Y(1:2),x3,Y(3:end),x8];
                                if not(sum(Y>1) + sum(Y<0))
                                    X{end+1} = Y;
                                end
                            end

                        else % This means NO(T1 fixed and S1 free).
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            Y = NT_Solver_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            % Y should have the values of T1, S1,
                            % T2, T3 and T4.  Y = [x2,x3,x4,x6,x7].
                            if isreal(Y)
                                x5 = (x8 - Y(3)*k(2) + k(1)) / k(3); % Y(3) = x4.
                                Y = [x1,Y(1:3),x5,Y(4:5),x8];
                                if not(sum(Y>1) + sum(Y<0))
                                    X{end+1} = Y;
                                end
                            end
                        end

                    elseif sets_aux(1) == 1 && bins_aux(1) == 1

                        % We ignore this case.

                    else % Here Phi_F is free.

                        if ~isempty(sets_aux) && sets_aux(1) == 2 && (length(sets_aux) == 1 || sets_aux(2) ~= 3) % This means T1 fixed and S1 free.
                            % Puedo ignorar c y nc.
                            nQ = nQ(2:end,2:end);
                            nb = nb(2:end);
                            nd = nd(2:end);
                            [~,nb,~,nd] = Mat_Red(nQ,nb,nc,nd,1,bins_aux(1)); % Reduzco T1.
                            T_Lambda_S1 = nb(1)/nd(1);
                            if (T_Lambda_F/T_Lambda_S1 > 0) && abs((T_Lambda_F/T_Lambda_S1) - 1) < Sens_Lambda
                                UseMatlabMin = 1;
                            end

                        else % NO(T1 fixed and S1 free).

                            nQ = nQ(2:end,2:end);
                            nb = nb(2:end);
                            nd = nd(2:end);
                            sets_aux = sets_aux - 1;
                            [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_F,sets_aux,bins_aux);
                            % Y should have the values of T1, S1,
                            % T3 and T4. Y = [x2,x3,x4,x6,x7].

                            if isreal(Y)
                                x5 = (x8 - Y(3)*k(2) + k(1)) / k(3); % Y(3) = x4.
                                Y = [0,Y(1:3),x5,Y(4:5)]; % Y = [0,x2,x3,x4,x5,x6,x7]. We add the 0 to do the next computation.
                                x1 = (Y*Q*Y'+Y*b+c)/(-b(1));
                                Y(1) = x1;
                                Y = [Y,x8];
                                if not(sum(Y>1) + sum(Y<0))
                                    X{end+1} = Y;
                                end
                            end

                        end

                    end % We finish checking Phi_F.

                else % 8 is not in sets(i,:) and {4,5} is not in sets(i,:).

                    if sets(i,1) == 1 && bins(j,1) == 0 % This means Phi_F = 0.
                        x1 = 0;
                        nQ = Q(2:end,2:end);
                        nb = b(2:end);
                        nc = c;
                        nd = d(2:end);
                        sets_aux = sets(i,2:end) - 1;
                        bins_aux = bins(j,2:end);

                        if ~isempty(sets_aux) && sets_aux(1) == 1 && (length(sets_aux) == 1 || sets_aux(2) ~= 2) % This means T1 fixed and S1 free.

                            if length(sets_aux) > 1 && sets_aux(2) == 3 && (length(sets_aux) == 2 || sets_aux(3) ~= 4) % This means T1 fixed and S1 free. % T1-F, S1-NF, T2-F and S2-NF.

                                [~,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                T_Lambda_S1 = nb(1)/nd(1);
                                T_Lambda_S2 = nb(2)/nd(2);
                                if ((T_Lambda_S1)/(T_Lambda_S2) > 0) &&  abs((T_Lambda_S1)/(T_Lambda_S2) - 1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end

                            else % T1 fixed, S1 free and NO(T2 fixed and S2 free).

                                x2 = bins_aux(1);
                                [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux(1),bins_aux(1)); % This reduces only T1.
                                sets_aux = sets_aux(2:end) - 1;
                                bins_aux = bins_aux(2:end);
                                T_Lambda_S1 = nb(1)/nd(1);
                                nb1 = nb(1);

                                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,1,0); % This removes S1.
                                sQ = nQ;
                                sb = nb;
                                sc = nc;
                                sets_aux = sets_aux - 1; % No S1 anymore.
                                [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_S1,sets_aux,bins_aux);
                                % Y should have the values of T2, S2,
                                % T3 and T4. Y = [x4,x5,x6,x7].

                                if isreal(Y)
                                    x3 = (Y*sQ*Y'+Y*sb+sc)/(-nb1);
                                    x8 = Y(2)*k(3) + Y(1)*k(2) - k(1); % Y(1) = x4, Y(2) = x5.
                                    Y = [x1,x2,x3,Y,x8];
                                    if not(sum(Y>1) + sum(Y<0))
                                        X{end+1} = Y;
                                    end
                                end

                            end

                        elseif ismember(3,sets_aux) && not(ismember(4,sets_aux))

                            [nQ,nb,nc,nd,sets_aux,bins_aux] = Change_Index(nQ,nb,nc,nd,sets_aux,bins_aux);

                            % Here we exchange the positions of
                            % T1<-->T2 and S1<-->S2 to make the system
                            % easier to solve.

                            x4 = bins_aux(1);
                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux(1),bins_aux(1)); % This reduces only T2.
                            sets_aux = sets_aux(2:end) - 1;
                            bins_aux = bins_aux(2:end);
                            T_Lambda_S2 = nb(1)/nd(1);
                            nb1 = nb(1);

                            [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,1,0); % This removes S2.
                            sQ = nQ;
                            sb = nb;
                            sc = nc;
                            sets_aux = sets_aux - 1; % No S2 anymore.

                            if length(nQ) ~= length(sets_aux)

                                [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                                Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_S2,sets_aux,bins_aux);
                                % Y should have the values of T1, S1,
                                % T3 and T4. Y = [x2,x3,x6,x7].

                                if isreal(Y)
                                    x5 = (Y*sQ*Y'+Y*sb+sc)/(-nb1); % This is S2.
                                    x8 = x5*k(3) + x4*k(2) - k(1);
                                    Y = [x1,Y(1:2),x4,x5,Y(3:4),x8];
                                    if not(sum(Y>1) + sum(Y<0))
                                        X{end+1} = Y;
                                    end
                                end

                            else % This case is m == 6 with sets(i,:) = [1,2,3,4,6,7].
                                % To avoid problems with the exchange,
                                % we directly compute all.

                                x2 = bins(j,sets(i,:)==2);
                                x3 = bins(j,sets(i,:)==3);
                                x6 = bins(j,sets(i,:)==6);
                                x7 = bins(j,sets(i,:)==7);
                                Y = [x1,x2,x3,x4,0,x6,x7];
                                x5 = (Y*Q*Y' + Y*b + c) / (2*Q(4,5)*x4);
                                x8 = x5*k(3) + x4*k(2) - k(1);
                                Y(5) = x5;
                                Y = [Y,x8];
                                if not(sum(Y>1) + sum(Y<0))
                                    X{end+1} = Y;
                                end

                            end

                        else % NO(T1 fixed and S1 free) and NO(T2 fixed and S2 free).

                            [nQ,nb,nc,nd] = Mat_Red_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            Y = NT_Solver_Bins(nQ,nb,nc,nd,sets_aux,bins_aux);
                            % Y should have T1, S1, T2, S2, T3 and T4.
                            % Y = [x2,x3,x4,x5,x6,x7].

                            if isreal(Y)
                                x8 = Y(4)*k(3) + Y(3)*k(2) - k(1); % Y(3) = x4 and Y(4) = x5.
                                Y = [x1,Y,x8];
                                if not(sum(Y>1) + sum(Y<0))
                                    X{end+1} = Y;
                                end
                            end
                        end

                    elseif sets(i,1) == 1 && bins(j,1) == 1

                        % We ignore this case.

                    else % Phi_F free.

                        sets_aux = sets(i,:);
                        bins_aux = bins(j,:);

                        if  ismember(2,sets_aux) && not(ismember(3,sets_aux))

                            if  ismember(4,sets_aux) && not(ismember(5,sets_aux))

                                [~,nb,~,nd] = Mat_Red_Bins(Q,b,c,d,sets_aux,bins_aux);

                                if (nb(1)/nd(1))/(nb(2)/nd(2)) > 0 && (nb(1)/nd(1))/(nb(3)/nd(3)) > 0 && abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) < ...
                                    Sens_Lambda && abs((nb(1)/nd(1))/(nb(3)/nd(3))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end

                            else % if not(ismember(4,sets_aux) && not(ismember(5,sets_aux)))

                                [~,nb,~,nd] = Mat_Red_Bins(Q,b,c,d,sets_aux,bins_aux);

                                if abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end

                            end

                        elseif ismember(4,sets_aux) && not(ismember(5,sets_aux)) % && not(ismember(2,sets_aux) && not(ismember(3,sets_aux)))

                            [~,nb,~,nd] = Mat_Red_Bins(Q,b,c,d,sets_aux,bins_aux);

                            if ismember(2,sets_aux) && ismember(3,sets_aux)

                                if abs((nb(1)/nd(1))/(nb(2)/nd(2))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end

                            elseif not(ismember(2,sets_aux)) && not(ismember(3,sets_aux))

                                if abs((nb(1)/nd(1))/(nb(4)/nd(4))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end

                            else % if not(ismember(2,sets_aux)) && ismember(3,sets_aux)

                                if abs((nb(1)/nd(1))/(nb(3)/nd(3))-1) < Sens_Lambda
                                    UseMatlabMin = 1;
                                end
                            end

                        else % Phi_F free and no other Lambdas.

                            [nQ,nb,nc,nd] = Mat_Red_Bins(Q,b,c,d,sets(i,:),bins(j,:));
                            [nQ,nb,~,nd] = Mat_Red_Bins(nQ,nb,nc,nd,1,0); % This removes Phi_F and reduces the system.
                            sets_aux = sets_aux - 1;

                            Y = T_Solver_Bins(nQ,nb,nd,T_Lambda_F,sets_aux,bins_aux);
                            % Y has T1, S1, T2, S2, T3 and T4. Y =
                            % [x2,x3,x4,x5,x6,x7].

                            if isreal(Y)
                                Y = [0,Y]; % We add this 0 only to do the next computation.
                                x1 = (Y*Q*Y'+Y*b+c)/(-b(1));
                                Y(1) = x1; % Y = [x1,x2,x3,x4,x5,x6,x7].
                                x8 = Y(5)*k(3) + Y(4)*k(2) - k(1);
                                Y = [Y,x8];
                                if not(sum(Y>1) + sum(Y<0))
                                    X{end+1} = Y;
                                end
                            end

                        end

                    end

                end

            end % For.
        end % For.

        if UseMatlabMin == 1

            options = optimoptions('fmincon','Display','notify');
            Minimum = Quad_FMC_Erik(Q,b,c,d,k,options)';
            Minimum(8) = Minimum(5)*k(3) + Minimum(4)*k(2) - k(1);
            disp('We had to use fmincon!');
            EndWhile = 1;

        elseif length(X) > 1

            Vals = zeros(1,length(X)-1);
            for i = 1:length(X)-1
                Vals(i) = X{i+1}(1:7)*d;
            end
            Minimum = X{find(Vals==min(Vals))+1};
            SolExists = 1;
            if WeIncreased == 1
                EndWhile = 1;
            end

        elseif length(X) == 1 && SolExists == 1

            EndWhile = 1;
            
        elseif length(X) == 1
            
            m = m + 2;
            WeIncreased = 1;

        end
            
    end % While.
    
%     fprintf('With m = %d and Demand = %d',m,-c);
    
end