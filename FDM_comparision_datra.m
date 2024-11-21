grid_sizes = [2, 4];
r = 16;
elapsed_times_warmup = zeros(length(grid_sizes), 1);
elapsed_times_post_warmup = zeros(length(grid_sizes), 1);
verification_errors = zeros(length(grid_sizes), 1);

for i = 1:length(grid_sizes)
    Nx = grid_sizes(i);
    Ny = grid_sizes(i)*r;
    Nz = grid_sizes(i)*r;
    
    Ax = tridiag(Nx); Ix = speye(Nx);
    Ay = tridiag(Ny); Iy = speye(Ny);
    Az = tridiag(Nz); Iz = speye(Nz);
    
    m = Nx * Ny * Nz;
    f = rand(m, 1);
    
    % Run FDM solver with warm-up
    [u, time_warmup] = FDM(Ax, Ay, Az, f);
    elapsed_times_warmup(i) = time_warmup;
    
    % Run FDM solver again after warm-up
    [u, time_post_warmup] = FDM(Ax, Ay, Az, f);
    elapsed_times_post_warmup(i) = time_post_warmup;
    
    % Verification
    A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);
    verification_errors(i) = norm(u - A_true \ f);
    
    % Display results for the current grid size
    fprintf('Grid size: %d x %d x %d\n', Nx, Ny, Nz);
    fprintf('Warm-up time: %.4f seconds\n', time_warmup);
    fprintf('Post-warm-up time: %.4f seconds\n', time_post_warmup);
    fprintf('Verification error: %.4e\n\n', verification_errors(i));
end

% Display summary table
fprintf('\nResults Summary:\n');
fprintf('Grid Size\tWarm-up Time (s)\tPost-Warm-up Time (s)\tVerification Error\n');
for i = 1:length(grid_sizes)
    fprintf('%dx%dx%d\t%.4f\t\t%.4f\t\t\t%.4e\n', ...
            grid_sizes(i), grid_sizes(i), grid_sizes(i), ...
            elapsed_times_warmup(i), elapsed_times_post_warmup(i), verification_errors(i));
end
