function rho = spectral_radius_var(A)

    [n, ~, p] = size(A);

    % reshape A into block row
    top = reshape(A, n, n*p);

    if p == 1
        C = top;
    else
        % companion matrix
        C = [top;
             eye(n*(p-1)), zeros(n*(p-1), n)];
    end

    % eigenvalues
    eigvals = eig(C);

    % spectral radius
    rho = max(abs(eigvals));
end