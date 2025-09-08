function f = chebfun3f_fixed_deg(fun, A, m)
    % fun: function handle, vectorized function of 3 variables
    % A: computational domain, matrix of size 3 x 2
    % m: degrees of Chebyshev polynomials for each dimension, vector of size 3
    % return a chebfun3 object
    f = chebfun3();
    
    A = parse_A(A,3);
    f.domain = [A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2)];
    % Set preferences:
    pref             = chebfunpref();
    tech             = pref.tech();
    prefStruct       = pref.cheb3Prefs;
    tpref            = tech.techPref;
    grid             = tpref.minSamples;
    maxSample        = tpref.maxLength;     % max polynomialDeg (not implemented)
    maxSamplePhase1  = 363;                 % max coarseResolution (not implemented)
    maxRank          = prefStruct.maxRank;  % max rank (not implemented)
    maxRestarts      = 10;                  % max restarts after the sample test
    pseudoLevel      = prefStruct.chebfun3eps;
    passSampleTest   = prefStruct.sampleTest;
    absTol           = pseudoLevel;


    chebX = @(i,m) A(1,1) + ((-cos((i-1).*pi/(m-1))) + 1)*(A(1,2)-A(1,1))/2;
    chebY = @(i,m) A(2,1) + ((-cos((i-1).*pi/(m-1))) + 1)*(A(2,2)-A(2,1))/2;
    chebZ = @(i,m) A(3,1) + ((-cos((i-1).*pi/(m-1))) + 1)*(A(3,2)-A(3,1))/2;

    n = 3;
    r = [6,6,6];
    restarts = 0; vectorize = 0;

    %% Phase 1: fiberwise ACA
    J = initializeIndexRandomly(r(2), m(2), restarts);
    K = initializeIndexRandomly(r(3), m(3), restarts);

    % Handle to evaluate tensor entries of T_c
    ff = @(i,j,k) fun(chebX(i,m(1)),chebY(j,m(2)),chebZ(k,m(3)));

    for iterations = 1:2
        % ACA 1
        T1 = evalTensor(1:m(1),J,K,ff,vectorize);

        % Does the function blow up or evaluate to NaN?:
        if ( isinf(max(abs(T1(:)))) )
            error('CHEBFUN:CHEBFUN3:chebfun3f:inf', ...
                'Function returned INF when evaluated');
        elseif ( any(isnan(T1(:))) )
            error('CHEBFUN:CHEBFUN3:chebfun3f:nan', ...
                'Function returned NaN when evaluated');
        end

        T1 = reshape(T1,m(1),r(2)*r(3));
        [~, absTol] = getTol(T1, pseudoLevel, absTol,A(1,2)-A(1,1),tech);
        [Uc, ~, ~, I,I2] = ACA(T1, absTol, m(1));
        r(1) = size(I,2);
        JT1 = J;
        KT1 = K;
        
        % ACA 2
        T2 = evalTensor(I,1:m(2),K,ff,vectorize);
        T2 = reshape(permute(T2,[2,1,3]),m(2),r(1)*r(3));
        [~, absTol] = getTol(T2, pseudoLevel, absTol,A(2,2)-A(2,1),tech);
        [Vc, ~, ~, J, J2] = ACA(T2, absTol, m(2));
        r(2) = size(J,2);
        KT2 = K;
        IT2 = I;
        
        % ACA 3
        T3 = evalTensor(I,J,1:m(3), ff,vectorize);
        T3 = reshape(permute(T3,[3,1,2]),m(3),r(1)*r(2));
        [relTol, absTol] = getTol(T3, pseudoLevel, absTol,A(3,2)-A(3,1),tech);
        [Wc, ~, ~, K, K2] = ACA(T3, absTol, m(3));
        r(3) = size(K,2);
        IT3 = I;
        JT3 = J;
    end

    if ( size(I,1) == 0 || size(J,1) == 0 || size(K,1) == 0 )
        f.cols = chebfun(zeros([n(1),1]), [A(1,1),A(1,2)], pref);
        f.rows = chebfun(zeros([n(2),1]), [A(2,1),A(2,2)], pref);
        f.tubes = chebfun(zeros([n(3),1]), [A(3,1),A(3,2)], pref);
        f.core = 0;
        return
    end

    %% Phase 2: Refinement (no refine here since m fixed)
    Uf = Uc;
    Vf = Vc;
    Wf = Wc;

    [~, absTol] = getTol(Uf, pseudoLevel, absTol, A(1,2)-A(1,1),tech);
    [~, absTol] = getTol(Vf, pseudoLevel, absTol, A(2,2)-A(2,1),tech);
    [~, absTol] = getTol(Wf, pseudoLevel, absTol, A(3,2)-A(3,1),tech);

    %% Phase 3: reconstruct the core tensor
    % Compute factor matrices
    [QU,~] = qr(Uf,0);
    [I, QUI] = DEIM(QU);
    [QV,~] = qr(Vf,0);
    [J, QVJ] = DEIM(QV);
    [QW,~] = qr(Wf,0);
    [K, QWK] = DEIM(QW);

    % Scaling to ensure the factor matrices contain decaying coefficients
    tmpCore  = invtprod(evalTensor(I,J,K,ff,vectorize), QUI,QVJ,QWK);
    colScaling = max(abs(tmpCore),[],[2,3]);
    rowScaling = max(abs(tmpCore),[],[1,3]);
    tubeScaling = squeeze(max(abs(tmpCore),[],[1,2]));
       
    % Store chebfun3 object
    f.cols  = chebfun(QU*diag(colScaling), [A(1,1),A(1,2)], pref);
    f.rows  = chebfun(QV*diag(rowScaling), [A(2,1),A(2,2)], pref);
    f.tubes = chebfun(QW*diag(tubeScaling), [A(3,1),A(3,2)], pref);
    f.core  = invtprod(tmpCore,diag(colScaling),diag(rowScaling),diag(tubeScaling));
end


function X = initializeIndexRandomly(r, maxVal, restarts)
    %  Random initialization of indices by drawing one index in each of r
    %  subintervals of equal length
    
    box = floor(maxVal/r);
    X = [];
    rngprev = rng();
    rng(16051821+restarts,'twister');             % pseudorandom
    for i = 1:r
        val = i*box + randi(box,1);
        X = [X,val];
    end
    rng(rngprev);
end

function T = evalTensor(I, J, K, ff,vectorize)
    % Evaluate the tensor ff at indices specified by I,J,K
    
    if ( vectorize == 0 ) % we can use the efficient evaluations
        n = [numel(I),numel(J),numel(K)];
        x = zeros([n(1),1,1]);
        x(:,1,1) = I;
        X = repmat(x,1,n(2),n(3));
        y = zeros([1,n(2),1]);
        y(1,:,1) = J;
        Y = repmat(y,n(1),1,n(3));
        z = zeros([1,1,n(3)]);
        z(1,1,:) = K;
        Z = repmat(z,n(1),n(2),1);
        if numel(X) > 0
            T = ff(X,Y,Z);
        else
            T = [];
        end
    else % we need for loops as f is not vectorizable
        T = zeros(size(I,2),size(J,2),size(K,2));
        for i = 1:size(I,2)
            for j =1:size(J,2)
                for k = 1:size(K,2)
                    T(i,j,k) = ff(I(i),J(j),K(k));
                end
            end
        end
    end
end

function [relTol, absTol] = getTol(M, pseudoLevel, tolOld,domDiff,tech)
    % Get suitable tolerances as in Chebfun3 (see
    % https://github.com/chebfun/chebfun/issues/1491)
    
    relTol = 2*size(M,1)^(4/5) * pseudoLevel;
    vscale = max(abs(M(:)));
    cheb = @(i,n) -cos((i-1).*pi/(n-1));
    if isa(tech,'trigtech')
        %use equispace points instead
        cheb = @(i,n) -1 +  (i-1)/(n) * 2;
    end
    points = 1:size(M,1);
    points = cheb(points, size(M,1));
    gradNorms = zeros([1,size(M,1)]);
    for i = 1:size(M,2)
        gradNorms(i) = max(abs(diff(M(:,i)) ./ diff(points)'));
    end
    gradNorms = max(gradNorms);
    absTol = max(max(domDiff.*gradNorms), vscale) * relTol;
    absTol = max([absTol, tolOld, pseudoLevel]);
end

function [Ac, At, Ar, rowInd, colInd] = ACA(A, tol, maxIter)
    % Adaptive Cross Approximation with full pivoting
    
    rowInd = [];
    colInd = [];
    Aoriginal = A;
    
    for iter = 1:maxIter
        
        [error,I2] = max(abs(A(:)));
        if ( isempty(error) || error < tol )
            Ac = Aoriginal(:,colInd);
            Ar = Aoriginal(rowInd,:)';
            At = Aoriginal(rowInd,colInd);
            return
        end
        
        [I,J] = ind2sub(size(A), I2);
        rowInd = [rowInd, I];
        colInd = [colInd, J];
        
        A = A-A(:,J)*A(I,:)./A(I,J);
    end
    Ac = Aoriginal(:,colInd);
    Ar = Aoriginal(rowInd,:)';
    At = Aoriginal(rowInd,colInd);
end

function [indices, UI] = DEIM(U)
    % Discrete Empirical Interpolation
    
    indices = [];
    [~, I] = max(abs(U(:,1)));
    indices = [indices,I];
    for l = 2:size(U,2)
        c = U(indices,1:(l-1)) \ U(indices,l);
        r = U(:,l) - U(:,1:(l-1))*c;
        [~, I] = max(abs(r));
        indices = [indices,I];
    end
    
    if ( nargout > 1 )
        UI = U(indices,:);
    end
    
end

function X = invtprod(X,U,V,W)
    % Evaluate X times_1 inv(U) times_2 inv(V) times_3 inv(W) using backslash
    
    n = [size(X,1),size(X,2),size(X,3)];
    m = [size(U,1),size(V,1),size(W,1)];
    X = reshape(U\reshape(X,[n(1),n(2)*n(3)]),[m(1),n(2),n(3)]);
    X = permute(reshape(V\reshape(permute(X,[2,1,3]),[n(2),m(1)*n(3)]),[m(2),m(1),n(3)]),[2,1,3]);
    X = permute(reshape(W\reshape(permute(X,[3,2,1]),[n(3),m(2)*m(1)]),[m(3),m(2),m(1)]),[3,2,1]);
    
end