function [X] = NT_Solver(Q,b,c,d)

    lambda = sqrt((b'*(Q\b)-4*c)/(d'*(Q\d)));
    if isreal(lambda)
        
        X = Q\(lambda*d-b)/2;
        A = null(X(:).');
        if all(eig(A'*Q*A) < 0)
            return;
        end
        
        lambda = -lambda;
        X = Q\(lambda*d-b)/2;
        A = null(X(:).');
        AllPos = all(eig(A'*Q*A) > 0);
        
        if AllPos
            return;
        end
        
        X = 1i;
        
    else
        X = 1i;
    end

end

% Here we check that lambda and the matrix A'QA have different signs,
% where for the matrix we use the positive or negative definition.