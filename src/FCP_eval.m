function vals = FCP_eval(lambda,Cheb,x,A)
    % evaluate the functional CP at x
    % lambda: weights of the FCP format, R x 1
    % Cheb: Chebyshev coefficients, cell of d matrices of size m_i x R
    % x: evaluation points, matrix of size N x d
    % A: computational domain, matrix of size d x 2
    % vals: values of the functional CP at x, N x 1
    [d,R,m] = parse_U(Cheb);
    A = parse_A(A,d);
    [x,N] = parse_x(x,d);

    % scale x to [-1,1]
    x = 2*(x-A(:,1)')./(A(:,2)'-A(:,1)')-1;
    
    % if m's are the same
    if all(m==m(1))
        % this is ~ 1 times faster than the following when m's are the same
        T = reshape((cos(acos(x(:))*(0:m(1)-1)))',[m(1),N,d]);
        T = pagetranspose(T);
        Cheb = reshape(cell2mat(Cheb),[m(1),R,d]);
        T_times_coeff = pagemtimes(T,Cheb);
        vals = prod(T_times_coeff,3) * lambda;
    else
        T = arrayfun(@(i,m_i) cos(acos(x(:,i))*(0:m_i-1)),1:d,m,'UniformOutput',false);
        T_times_coeff = cellfun(@(T_i,Cheb_i) T_i*Cheb_i,T,Cheb,'UniformOutput',false);
        vals = prod(reshape(cell2mat(T_times_coeff),[N,R,d]),3) * lambda;
    end
end

function [x,N] = parse_x(x,d)
    % check the input is of the correct format
    % if x is a vector of length d, then convert it to a matrix of size 1 x d
    if isvector(x) 
        if length(x)~=d
            error('The evaluation points should be a matrix of size N x d.');
        end
        x = x(:)';
        N = 1;
    elseif ismatrix(x)
        if size(x,2)~=d
            error('The evaluation points should be a matrix of size N x d.');
        end
        N = size(x,1);
    else
        error('The evaluation points should be a matrix of size N x d.');
    end
end