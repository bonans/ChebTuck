function X_sub = CP_get_subtensor(lambda, U, idx)
    % Get the subtensor in full format of a tensor given by its CP decomposition
    % lambda: weights of the rank-1 tensors, R x 1
    % U: factor matrices, cell of d matrices of size n_i x R, i.e., 
    % the original tensor is of size n_1 x n_2 x ... x n_d
    % idx: indices of the subtensor, cell of d vectors of size m_i
    % X_sub: subtensor of size m_1 x m_2 x ... x m_d

    % get the subtensor in CP format
    X_sub_CP = cellfun(@(U, idx) U(idx, :), U, idx, 'UniformOutput', false);

    if TensorToolboxInstalled()
        X_sub = full(ktensor(lambda, X_sub_CP)).data;
    else
        % otherwise, compute the subtensor in full format
        X_sub = lambda(1) * rank1_2_full(cellfun(@(x) x(:, 1), X_sub_CP, 'UniformOutput', false));
        for r = 2:length(lambda)
            X_sub = X_sub + lambda(r) * rank1_2_full(cellfun(@(x) x(:, r), X_sub_CP, 'UniformOutput', false));
        end
    end
end