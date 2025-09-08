function [Cheb,cheb_errors] = CP_2_FCP(grids,CU,A,m)
% convert CP to functional CP using Chebyshev interpolation
% grids: grid points, cell of d vectors of size n_i x 1
% CU: factor matrices, cell of d matrices of size n_i x R
% A: computational domain, matrix of size d x 2
% m: degrees of Chebyshev polynomials for each dimension, vector of size d


% process input
[d,R,n] = parse_U(CU);
grids = parse_grids(grids,d,n);
A = parse_A(A,d);
m = parse_m(m,d);

Cheb = cell(1,d);
cheb_errors = zeros(d,R);
for ith = 1:d
    f_cheb = chebfun(@(x) interp1(grids{ith}',CU{ith},x,'spline','extrap'),A(ith,:),m(ith));
    cheb_errors(ith,:) = vecnorm(f_cheb(grids{ith})-CU{ith},'inf');
    Cheb{ith} = chebcoeffs(f_cheb);
end

end

function m = parse_m(m,d)
    % check the input is of the correct format
    % if it is a number, assume the degree is the same for each dimension
    if isscalar(m)
        m = repmat(m,d,1);
    elseif isvector(m) && length(m)~=d
       error('The degrees of Chebyshev polynomials m should be a vector of size d.'); 
    end
end