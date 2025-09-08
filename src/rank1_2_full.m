function X_full = rank1_2_full(U)
    % Get the full tensor from rank-1 tensors
    % U: factor matrices, cell of d matrices of size n_i x 1, i.e.,
    % the original tensor is of size n_1 x n_2 x ... x n_d
    % X_full: full tensor of size n_1 x n_2 x ... x n_d
    
    d = length(U);
    n = cellfun(@(x) size(x, 1), U);
    n = reshape(n,1,[]);
    
    X_full = U{d};
    for ii = d-1:-1:1
        X_full = kron(X_full, U{ii});
    end
    X_full = reshape(X_full, n);
end