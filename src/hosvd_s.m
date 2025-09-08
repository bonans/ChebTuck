function [A, U] = hosvd_s(A, tol, varargin)
    % compute the HOSVD of a tensor with a specified tolerance or rank
    % A: tensor of size n1 x n2 x ... x nd
    % tol: tolerance for the approximation
    % varargin: ranks of the approximation, vector of size d

    % process input
    d = ndims(A);

    if nargin == 2
        r = ones(1,d);
        given_r = false;
        reltol = tol * norm(A,'fro')/sqrt(d);
    elseif nargin == 3
        r = varargin{1};
        if length(r) ~= d
            error('Number of ranks should equal the dimension of A.');
        end
        given_r = true;
    else
        error('The number of input arguments should be 2 or 3.');
    end

    U = cell(d,1);
    % compute the factor matrices
    for k = 1:d
        Ak = ten2mat(A,k);
        if given_r
            [U{k},~,~] = svds(Ak,r(k));
            
        else
            [Uk,Sk,~] = svd(Ak,"econ");
            Sk = diag(Sk);
            Sksum = sqrt(cumsum(Sk.^2, "reverse"));
            r(k) = find(Sksum > reltol, 1, 'last');
            U{k} = Uk(:,1:r(k));
        end
    end

    % compute the core tensor
    A = tucker2full(A, cellfun(@(x) x', U, 'UniformOutput', false));
end