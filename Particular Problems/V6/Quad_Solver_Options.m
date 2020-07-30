function [X] = Quad_Solver_Options(iQ,ib,ic,id,ik,options,H)



    lb = zeros(1,length(iQ));
    ub = ones(1,length(iQ));

    fun = @(x)quadobj(x,0,id,0);
    nonlconstr = @(x)quadconstr(x,H,k,d);
    x0 = zeros(length(iQ),1); % Column vector.
%     [X,fval,eflag,output,lambda] = fmincon(fun,x0,...
%         [],[],[],[],lb,ub,nonlconstr,options);
    [X,~,~,~,~] = fmincon(fun,x0,...
        [],[],[],[],lb,ub,nonlconstr,options);
        
end
