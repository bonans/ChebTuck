function P = uni2cheb(dom, n, m)
    % the matrix P of size m x n such that
    % Pa, where a is a vector of length n containing the values at 
    % the uniform grid points, is the extrapolated values
    % at the Chebyshev grid, by using the spline interpolation

    % dom = [a, b] is the domain of the uniform grid
    % n is the number of uniform grid points
    % m is the number of Chebyshev grid points

    unigrid = linspace(dom(1), dom(2), n+1);
    unigrid = unigrid(1:end-1); unigrid = unigrid(:);
    chebgrid = chebpts(m,[dom(1), dom(2)]); chebgrid = chebgrid(:);
    % P = zeros(m, n);
    % for ii = 1:n
    %     vals = zeros(n, 1);
    %     vals(ii) = 1;
    %     P(:, ii) = interp1(unigrid,vals,chebgrid,'spline','extrap');
    % end
    P = interp1(unigrid,eye(n),chebgrid,'linear','extrap');
end