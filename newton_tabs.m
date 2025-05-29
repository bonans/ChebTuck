close all;
d = 3; A = 1;
% create the table of approximation errors
tol = 1e-7;
nlist = [256, 512, 1024, 2048, 4096];
num_n = length(nlist);
mlist = [129,257,513,1025,2049,4097,8193,16385];
num_m = length(mlist);
errors = zeros(num_n,num_m);

errors_long = zeros(num_n,num_m);
TuckerRanks_long = zeros(num_n,num_m);

for ii = 1:num_n
    for jj = 1:num_m
        n1 = nlist(ii);
        m = mlist(jj) * ones(1,d);
        h1 = 2/n1;
        xcol = -A +0.5*h1:h1:A -0.5*h1;
        filename = ['data/Newton_CP_-1_1/Newt_canon_' num2str(n1) '.mat'];
        load(filename)

        % get the CP tensor
        U1r = U1r./h1;
        R = size(U1r,2);
        xi = ones(R,1);
        U = {U1r, U1r, U1r};

        % compute the middle slice of the tensor
        ns1 = floor(n1/2);
        F = CP_get_subtensor(xi,U,{1:n1,1:n1,ns1:ns1+1});
        X_slice = F(:,:,1);

        % compute ChebTuck
        f = ChebTuck({xi,U},m,[],tol);

        % compute the ChebTuck error
        [xx, yy, zz] = ndgrid(xcol,xcol,xcol(ns1:ns1+1));
        X_cheb1 = f(xx,yy,zz);
        X_ChebTuck = X_cheb1(:,:,1);

        errors(ii,jj) = norm(X_slice(:) - X_ChebTuck(:),inf)/norm(X_slice(:),inf);
        
        % long-range part
        Rl = ceil(R/2)-3;
        xi = ones(Rl,1);
        U = {U1r(:,1:Rl), U1r(:,1:Rl), U1r(:,1:Rl)};

        % compute the middle slice of the tensor
        ns1 = floor(n1/2);
        F = CP_get_subtensor(xi,U,{1:n1,1:n1,ns1:ns1+1});
        X_slice = F(:,:,1);

        % compute ChebTuck
        f = ChebTuck({xi,U},m,[],tol);

        % compute the ChebTuck error
        [xx, yy, zz] = ndgrid(xcol,xcol,xcol(ns1:ns1+1));
        X_cheb1 = f(xx,yy,zz);
        X_ChebTuck = X_cheb1(:,:,1);

        errors_long(ii,jj) = norm(X_slice(:) - X_ChebTuck(:),inf)/norm(X_slice(:),inf);
        TuckerRanks_long(ii,jj) = max(size(f.core));
    end
end