function A = parse_A(A,d)
    % parse the input for the computational domain
    % if it is a number, assume the domain is [-A,A] for each dimension
    % if it is a 2-vector, assume the domain is [A(1),A(2)] for each dimension
    if isscalar(A)
        A = repmat([-A,A],d,1);
    elseif isvector(A) && length(A)==2
        A = repmat([A(1),A(2)],d,1);
    elseif ismatrix(A) && any(size(A)~=[d,2])
        error('The computational domain should be a matrix of size 2 x d.');
    end
end