function vals = fun_grid2val(xx,yy,zz,F,grids)
    % given discretized function values on a grid, evaluate the function at a point in the domain
    % xx,yy,zz: points to evaluate the function, should be of the same size
    % F: function values, array of size n1 x n2 x n3
    % grids: grid points, cell of 3 vectors of size n_i x 1
    % A: computational domain, matrix of size 3 x 2

    %vals = interp3(grids{1},grids{2},grids{3},F,xx,yy,zz,'nearest');
    vals = interp3(grids{1},grids{2},grids{3},F,xx,yy,zz,'splines');
    % global num_evals;
    % num_evals = num_evals + numel(vals);
end