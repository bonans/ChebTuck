function X = ten_x_mat(X, U, inds)
    % tensor times matrices
    % X: tensor of size n1 x n2 x ... x nd
    % U: cell of length s matrices of size m_{i_k} x n_{i_k}
    % inds: vector of length s, the modes to multiply
    % return a tensor multiplied by U's at the specified modes in inds

    s = length(inds);
    d = ndims(X);
    if s == 1
        X = permute(tensorprod(X, U, inds, 2), [1:inds-1, d, inds:d-1]);
    elseif s > 1
        for i = 1:length(inds)
            X = permute(tensorprod(X, U{i}, inds(i), 2), [1:inds(i)-1, d, inds(i):d-1]);
        end
    end
end