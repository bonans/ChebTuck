function Xk = ten2mat(X, k)
    % tensor mode-k unfolding
    % X: array of size n1 x n2 x ... x nd
    % k: integer, mode-k unfolding
    % return a matrix of size nk x (n1*...*nk-1*nk+1*...*nd)
    d = ndims(X);
    n = size(X);
    X = permute(X, [k, 1:k-1, k+1:d]);
    Xk = reshape(X, n(k), []);
end