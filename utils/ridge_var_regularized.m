function [A, SIG] = ridge_var_regularized(X, p, lambda)

[n, T] = size(X);

% Build regression matrices
Y = X(:, p+1:T);
Z = [];

for k = 1:p
    Z = [Z; X(:, p+1-k:T-k)];
end

% Ridge solution
I = eye(size(Z,1));

A_big = (Y * Z') / (Z * Z' + lambda * I);

% Reshape into VAR form
A = zeros(n, n, p);
for k = 1:p
    A(:,:,k) = A_big(:, (k-1)*n+1:k*n);
end

% Residuals
E = Y - A_big * Z;
SIG = (E * E') / size(E,2);

end