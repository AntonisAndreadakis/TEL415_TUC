function X = real_colored_noise(m, C, N, method)
    if method == 1      % eigen
        [Q, L] = eig(C);
        F = Q * L^(1/2);
        res = F;
    elseif method ==2       % cholensky
        L = chol(C);
        res = L;
    else
        fprintf('Invalid method. Choose 1 for eig or 2 for chol.\n');
    end

    n = length(m);
    Y = randn(n, N);
    Z = res * Y;
    X = Z + m * ones(1,N);
end