function [f,ft] = chebfun3h_fixed_deg(fun, A, m, tol, varargin)
    % construct Hybrid Cheb-Tuck apprximation by 
    % 1. Construct the full 3D Chebyshev coefficient tensor C
    % 2. Apply HOSVD to the tensor C
    % 3. Return the approximated function in a chebfun3 object
    % 
    % fun: function handle, vectorized function of 3 variables
    % A: computational domain, matrix of size 3 x 2
    % m: degrees of Chebyshev polynomials for each dimension, vector of size 3
    % tol: tolerance of
    % return a chebfun3 object

    f = chebfun3();
    
    A = parse_A(A,3);   
    f.domain = [A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2)];

    xx = chebpts(m(1),[A(1,1),A(1,2)]);
    yy = chebpts(m(2),[A(2,1),A(2,2)]);
    zz = chebpts(m(3),[A(3,1),A(3,2)]);
    [xx,yy,zz] = ndgrid(xx,yy,zz);

    Fvals = fun(xx,yy,zz);

    Fcoefs = chebfun3t.vals2coeffs(Fvals);

    ft = chebfun3t();
    ft.domain = f.domain;
    ft.coeffs = Fcoefs;
    
    if isempty(varargin)
        [core, U] = hosvd_s(Fcoefs, tol);
    elseif length(varargin) == 1
        [core, U] = hosvd_s(Fcoefs, tol, varargin{1});
    else
        error('The number of input arguments should be 4 or 5.');
    end

    
    f.cols = chebfun(U{1}, A(1,:), 'coeffs');
    f.rows = chebfun(U{2}, A(2,:), 'coeffs');
    f.tubes = chebfun(U{3}, A(3,:), 'coeffs');
    f.core = core;

end