function [d,R,n] = parse_U(CU)
    % check the input for the factor matrices
    d = length(CU);
    Rs = cellfun(@(x) size(x,2),CU);
    if length(unique(Rs))>1
        error('The second dimension of each factor matrix should be the same.');
    end
    R = Rs(1);
    n = cellfun(@(x) size(x,1),CU);
end