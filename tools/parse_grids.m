function grids = parse_grids(grids,d,n)
    % check the input is of the correct format
    % if grids is not a cell, then convert it to a cell, each cell is the same vector
    if ~iscell(grids)
        if length(unique(n))>1
            error('The grid points should be a cell of d vectors of size n_i x 1.');
        end
        grids = repmat({grids},1,d);
    else
        if length(grids)~=d
            error('The grid points should be a cell of d vectors of size n_i x 1.');
        end
        if any(cellfun(@(x) length(x),grids)~=n)
            error('The grid points should be a cell of d vectors of size n_i x 1.');
        end
    end
    % check if the grid points are a column vector, otherwise convert it to a column vector
    if any(cellfun(@(x) size(x,2)~=1,grids))
        grids = cellfun(@(x) x(:),grids,'UniformOutput',false);
    end
end