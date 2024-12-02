% After assembling K and F in previous code
% Solve using direct method
d_temp = K \ F;
disp_direct = [d_temp; g]; % Append Dirichlet BC

% Solve using GMRES with different tolerances
tol_values = [1e-2, 1e-4, 1e-6];
maxit = 10000;
restart = [];

for tol = tol_values
    [d_temp_gmres, flag, relres, iter] = gmres(K, F, restart, tol, maxit);
    disp_gmres = [d_temp_gmres; g]; % Append Dirichlet BC

    % Compute the relative error between GMRES and direct solution
    rel_error = norm(disp_gmres - disp_direct) / norm(disp_direct);
    fprintf('GMRES with tol=%e: relative error = %e, iterations = %d\n', tol, rel_error, sum(iter));

    % You can also compute errors as before if needed
end
