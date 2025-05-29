function CP = CP_normalize(CP)
    % normalize the CP tensor
    % CP_in should be {lambda,U}
    % lambda is a vector of (positive) weights, R x 1
    % U is a cell of factor matrices, i.e. d matrices of size n_i x R
    % return a normalized CP tensor, i.e., 
    % the weights should be positive and 
    % the each column of each factor matrix should have unit norm

    lambda = CP{1};
    U = CP{2};
    [d,R] = parse_U(U); 

    % first multiply the weights back to U{1}
    U{1} = U{1} .* reshape(lambda,1,R);
    % then extract the norm of each column
    lambda = ones(1,R);
    for i = 1:d
        vecnormi = vecnorm(U{i});
        lambda = lambda .* vecnormi;
        U{i} = U{i} ./ vecnormi;
    end
    CP = {lambda',U};
end
