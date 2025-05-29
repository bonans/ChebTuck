function X = mat2ten(Xk, n, k)
    % given a matrix Xk of size nk x (n1*...*nk-1*nk+1*...*nd)
    % return a tensor of size n1 x n2 x ... x nd
    d = length(n);
    X = reshape(Xk, [n(k), n(1:k-1), n(k+1:d)]);
    X = ipermute(X, [k, 1:k-1, k+1:d]);
end