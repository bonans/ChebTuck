function vals = fun_gridCP2val(xx,yy,zz,lambda,U,grids,A, varargin)
    % given discretized function values on a grid in CP format, evaluate the
    % function at another grid
    % xx,yy,zz: the tensor product grid points to be evaluated, vectors of size m_i,
    % i.e., the returned vals is of size m_1 x m_2 x m_3
    % lambda: weights of the rank-1 tensors, R x 1
    % U: factor matrices, cell of 3 matrices of size n_i x R, i.e., 
    % the original tensor is of size n_1 x n_2 x n_3
    % grids: grid points, cell of 3 vectors of size n_i x 1
    % A: computational domain, matrix of size 3 x 2
    % Optional input:
    % interpolation type, default is 'spline'

    % parse optional inputs
    if isempty(varargin)
        interptype = 'spline';
    else
        interptype = varargin{1};
    end
    %interptype = 'nearest'; 
    [d,R,n] = parse_U(U);
    grids = parse_grids(grids,d,n);
    A = parse_A(A,d);

    % if xx,yy,zz are full grid, retrieve the grid vectors first
    if ndims(xx) == 3
        xx = reshape(xx(:,1,1),[],1);
        yy = reshape(yy(1,:,1),[],1);
        zz = reshape(zz(1,1,:),[],1);
    end


    
    % f_chebx = chebfun(@(x) interp1(grids{1}',U{1},x,interptype,'extrap'),A(1,:),m(1));
    % f_cheby = chebfun(@(x) interp1(grids{2}',U{2},x,interptype,'extrap'),A(2,:),m(2));
    % f_chebz = chebfun(@(x) interp1(grids{3}',U{3},x,interptype,'extrap'),A(3,:),m(3));

    CU = {interp1(grids{1}',U{1},xx,interptype,'extrap'), ...
            interp1(grids{2}',U{2},yy,interptype,'extrap'), ...
            interp1(grids{3}',U{3},zz,interptype,'extrap')};
    vals = CP_get_subtensor(lambda, CU, {1:length(xx), 1:length(yy), 1:length(zz)});
end