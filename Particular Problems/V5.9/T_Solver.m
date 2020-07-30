function [X] = T_Solver(Q,b,d,T_lambda)

    X = Q\(T_lambda*d-b)/2;
    A = null(X(:).');
    if (all(eig(A'*Q*A) < 0) && T_lambda > 0) || (all(eig(A'*Q*A) > 0) && T_lambda < 0)
        return;
    else
        X = 1i;
    end

end