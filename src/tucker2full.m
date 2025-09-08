function X = tucker2full(X, U)
    % compute the full tensor from the Tucker representation
    % X: core tensor of size r1 x r2 x ... x rd
    % U: cell array of length d of factor matrices of size ni x ri
    % return X: full tensor of size n1 x n2 x ... x nd
    
    for k = 1:length(U)
        X = tensorprod(X, U{k}, 1, 2);
    end

end