function M = Rand_Sym(n)

   d = 1*rand(n,1); % The diagonal values
   t = triu(bsxfun(@min,d,d.').*rand(n),1); % The upper trianglar random values
   M = diag(d)+t+t.'; % Put them together in a symmetric matrix
   
end