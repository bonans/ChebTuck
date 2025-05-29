function [f,ff] = ChebTuck(fun,m,A,tol,r,type)
    % construct ChebTuck approximation of a function/tensor
    % 3 cases: 
    % 1. fun is a function handle
    % 2. fun is a full tensor
    % 3. fun is a CP tensor
    % 
    % case 1: fun is a function handle
    % f = ChebTuck(fun) returns the ChebTuck approximation of fun.
    % fun is a (vectorized) function handle of three variables.
    % f is a chebfun3 object.
    %
    % f = ChebTuck(fun,m,A,tol) additionally specifies the Chebyshev
    % degrees m, the computational domain A, the HOSVD tolerance tol
    % the default values are
    % m = [100,100,100], A = [-1,1;-1,1;-1,1], tol = 1e-7
    %
    % f = ChebTuck(fun,m,A,tol,r) specifies the Tucker ranks r
    % this will override the tolerance tol. One can also use
    % f = ChebTuck(fun,[],[],[],r) to specify the ranks only.
    %
    % [f,ff] = ChebTuck(___) also returns the full Chebyshev interpolant
    % ff without Tucker approximation. ft is a chebfun3t object.
    %
    %
    % case 2: fun is a full tensor of size n1 x n2 x n3
    % the usage is the same as case 1 but the function handle 
    % is replaced by the full tensor. Note that this is very
    % expensive and only tractable for small sizes.
    %
    % case 3: fun is a CP tensor, i.e., {lambda,U}.
    % lambda is a vector of (positive) weights, R x 1
    % U is a cell of factor matrices, i.e. d matrices of size n_i x R
    %
    % [f,ff] = ChebTuck(fun,m,A,tol,r) usage is the same as case 1 but
    % tol is the RHOSVD tolerance for the Tucker approximation 
    % and ff is a CP tensor of the Chebyshev coefficients, i.e.,
    % a cell of factor matrices of size m_i x R
    % 
    % [_] = ChebTuck(___,type) additionally specifies the type of the approximation.
    % This is only relevant for case 3. type = 'R' (default) or 'F'.
    % 'F' computes the Chebyshev coefficients in full format and then 
    % use the HOSVD. 'R' directly computes the Tucker approximation from
    % CP tensor using the RHOSVD.
    % 
    % case 1, 2 corresponds to Alg. 1,2 in the paper.
    % case 3 with 'F', 'R' (default) corresponds to Alg. 3,4 in the paper.
    %
    
    % parse input
    if nargin < 2 || isempty(m)
        m = [100,100,100];
    end
    if nargin < 3 || isempty(A)
        A = [-1,1;-1,1;-1,1];
    end
    if nargin < 4 || isempty(tol)
        tol = 1e-7;
    end
    given_r = false;
    if nargin > 4 && ~isempty(r)
        given_r = true;
    end
    if nargin < 6 || isempty(type)
        type = 'R';
    end

    f = chebfun3();
    f.domain = [A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2)];

    if isa(fun,'float')
        n = size(fun); d = length(n);
    elseif iscell(fun)
        lambda = fun{1}; U = fun{2}; n = cellfun(@(x) size(x,1),U);
        [d,R,n] = parse_U(U); 
    end
    % uniform grid in the computational domain
    xx = (-1+1/n(1)):2/n(1):(1-1/n(1));
    yy = (-1+1/n(2)):2/n(2):(1-1/n(2));
    zz = (-1+1/n(3)):2/n(3):(1-1/n(3));
    grids = {xx,yy,zz};

    % case 1, 2, 3 with 'F': convert fun to a function handle
    if isa(fun,'float')
        % case 2: fun is a full tensor, Algorithm 2
        fun = @(x,y,z) interp3(xx,yy,zz,fun,x,y,z,'splines');
    elseif iscell(fun) && strcmp(type,'F')
        % case 3: fun is a CP tensor HOSVD, Algorithm 3
        fun = @(x,y,z) fun_gridCP2val(x,y,z,lambda,U,grids,A);
    end


    if isa(fun,'function_handle')
        % case 1: fun is a function handle

        % compute the function evaluations T
        xx = chebpts(m(1),[A(1,1),A(1,2)]);
        yy = chebpts(m(2),[A(2,1),A(2,2)]);
        zz = chebpts(m(3),[A(3,1),A(3,2)]);
        [xx,yy,zz] = ndgrid(xx,yy,zz);
        T = fun(xx,yy,zz);
        % compute the Chebyshev coefficients C
        C = chebfun3t.vals2coeffs(T);

        % return the full Chebyshev interpolant
        ff = chebfun3t();
        ff.domain = f.domain;
        ff.coeffs = C;

        % compute the Tucker decomposition of C
        if given_r
            [core,V] = hosvd_s(C,tol,r);
        else
            [core,V] = hosvd_s(C,tol);
        end
    elseif iscell(fun) && strcmp(type,'R')
        % case 3: fun is a CP tensor RHOSVD, Algorithm 4

        % compute the CP format of the Chebyshev coefficients
        
        CU = {chebtech2.vals2coeffs(interp1(grids{1}',U{1},chebpts(m(1),[A(1,1),A(1,2)]),'spline','extrap')), ...
            chebtech2.vals2coeffs(interp1(grids{2}',U{2},chebpts(m(2),[A(2,1),A(2,2)]),'spline','extrap')), ...
            chebtech2.vals2coeffs(interp1(grids{3}',U{3},chebpts(m(3),[A(3,1),A(3,2)]),'spline','extrap'))};
        
        % CU = cell(1,d);
        % for ith = 1:d
        %     f_cheb = chebfun(@(x) interp1(grids{ith}',U{ith},x,'spline','extrap'),A(ith,:),m(ith));
        %     CU{ith} = chebcoeffs(f_cheb);
        % end

        % return the CP tensor of the Chebyshev coefficients
        ff = CU;

        % compute the RHOSVD of the CP tensor
        if given_r
            [core,V] = rhosvd(lambda,CU,tol,r);
        else
            [core,V] = rhosvd(lambda,CU,tol);
        end
    end
    % return the ChebTuck format
    f.cols = chebfun(V{1}, A(1,:), 'coeffs');
    f.rows = chebfun(V{2}, A(2,:), 'coeffs');
    f.tubes = chebfun(V{3}, A(3,:), 'coeffs');
    f.core = core;
end