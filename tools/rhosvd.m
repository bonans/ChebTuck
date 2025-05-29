function [A, V] = rhosvd(lambda, U, tol, varargin)
    % compute the reduced HOSVD of a CP tensor with a specified tolerance or rank
    % lambda: weights of the rank-1 tensors, R x 1
    % U: factor matrices, cell of d matrices of size n_i x R, i.e., 
    % the original tensor is of size n_1 x n_2 x ... x n_d
    % tol: tolerance for the approximation
    % varargin: ranks of the approximation, vector of size d

    % process input
    [d,R,n] = parse_U(U);

    if nargin == 3
        r = ones(1,d);
        given_r = false;
        reltol = tol/d;
    elseif nargin == 4
        r = varargin{1};
        if length(r) ~= d
            error('Number of ranks should equal the dimension of A.');
        end
        given_r = true;
    else
        error('The number of input arguments should be 2 or 3.');
    end

    V = cell(d,1);
    Z = cell(d,1);
    % compute the factor matrices
    for k = 1:d
        if given_r
            [V{k},UkS,UkV] = svds(U{k}, r(k));
            Z{k} = UkS*UkV';
        else
            [UkU,UkS,UkV] = svd(U{k},"econ");
            Sk = diag(UkS);
            Sksum = sqrt(cumsum(Sk.^2, "reverse"));
            r(k) = find(Sksum > reltol, 1, 'last');
            V{k} = UkU(:,1:r(k));
            Z{k} = UkS(1:r(k),1:r(k))*UkV(:,1:r(k))';
        end
    end

    % compute the core tensor
    %lambda_diag = zeros(R,R,R);
    %lambda_diag(1:R^2+R+1:end) = lambda;
    %A = tucker2full(lambda_diag,Z);
    % lambda_diag = sptendiag(lambda,[R,R,R]);
    % A = ttm(lambda_diag,Z,1:d).data;
    A = CP_get_subtensor(lambda, Z, {1:r(1); 1:r(2); 1:r(3)});
end