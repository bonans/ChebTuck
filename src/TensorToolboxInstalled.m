function v = TensorToolboxInstalled()
    % Check if the Tensor Toolbox is installed
    p = path();
    v1 = contains(p,'tensor_toolbox');
    test_funs = {'cp_arls','symktensor','eig_sshopm','sptensor','ttensor','ktensor','tensor'};
    v2 = all(cellfun(@(x) exist(x,'file'),test_funs) > 0);
    v = v1 || v2;
end