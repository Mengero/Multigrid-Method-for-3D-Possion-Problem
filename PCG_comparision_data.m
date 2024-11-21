grid_sizes = [4, 8, 16];
r = 16;
elapsed_times_pcg = zeros(length(grid_sizes), 1);
iterations_pcg = zeros(length(grid_sizes), 1);
residual_norms_pcg = zeros(length(grid_sizes), 1);

elapsed_times_cg = zeros(length(grid_sizes), 1);
iterations_cg = zeros(length(grid_sizes), 1);
residual_norms_cg = zeros(length(grid_sizes), 1);

for i = 1:length(grid_sizes)
    Nx = grid_sizes(i);
    Ny = grid_sizes(i)*r;
    Nz = grid_sizes(i)*r;
    
    Ax = tridiag(Nx); Ix = speye(Nx);
    Ay = tridiag(Ny); Iy = speye(Ny);
    Az = tridiag(Nz); Iz = speye(Nz);
    
    m = Nx * Ny * Nz;
    f = rand(m, 1);
    u = zeros(size(f));
    A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);
    
    mode = 2;
    nsmooth = 50;

    % Preconditioned Conjugate Gradient (PCG) solver
    tic;
    M = @(x) multigridPreconditioner(A_true, Nx, Ny, Nz, x, mode, nsmooth);
    [x, iter, resl] = myPCG(A_true, f, 1e-8, 2000, M);
    elapsed_times_pcg(i) = toc;
    iterations_pcg(i) = iter;
    residual_norms_pcg(i) = norm(A_true * x - f);
    
    % Display results for the current grid size
    fprintf('Grid size: %d x %d x %d\n', Nx, Ny, Nz);
    fprintf('PCG: Elapsed time: %.4f seconds, Iterations: %d, Final residual norm: %.4e\n', ...
            elapsed_times_pcg(i), iterations_pcg(i), residual_norms_pcg(i));
end

% Display summary table
fprintf('\nResults Summary:\n');
fprintf('Grid Size\tPCG Time (s)\tPCG Iterations\tPCG Residual\tCG Time (s)\tCG Iterations\tCG Residual\n');
for i = 1:length(grid_sizes)
    fprintf('%dx%dx%d\t%.4f\t\t%d\t\t%.4e\t%.4f\t\t%d\t\t%.4e\n', ...
            grid_sizes(i), grid_sizes(i), grid_sizes(i), ...
            elapsed_times_pcg(i), iterations_pcg(i), residual_norms_pcg(i), ...
            elapsed_times_cg(i), iterations_cg(i), residual_norms_cg(i));
end
