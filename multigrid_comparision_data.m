grid_sizes = [2, 4, 16];
r = 16;
elapsed_times = zeros(length(grid_sizes), 1);
residual_norms = zeros(length(grid_sizes), 1);
loop_counts = zeros(length(grid_sizes), 1); % To store the number of while loop iterations

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
    
    residue = f;
    u = zeros(size(f));
    
    mode = 2;
    nsmooth = 50;
    
    beta_list = beta_list_generation(nsmooth);
    
    loop_count = 0; % Counter for the while loop
    
    tic;
    while norm(residue)/norm(f) > 1e-8
        [u, rr] = vv_cycle(A_true, Nx, Ny, Nz, f, u, mode, nsmooth, beta_list);
        residue = rr;
        loop_count = loop_count + 1; % Increment the loop counter
        disp(['Grid size ', num2str(grid_sizes(i)), ': Residual norm = ', num2str(norm(residue))]);
    end
    elapsed_time = toc;
    
    elapsed_times(i) = elapsed_time;
    residual_norms(i) = norm(A_true * u - f);
    loop_counts(i) = loop_count; % Store the number of while loop iterations
    
    disp(['Grid size ', num2str(grid_sizes(i)), '^3:']);
    disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
    disp(['Final residual norm: ', num2str(residual_norms(i))]);
    disp(['Number of while loop iterations: ', num2str(loop_count)]);
end

% Display final results
disp('Results:');
disp('Grid Size   Elapsed Time (s)   Final Residual Norm   Loop Count');
for i = 1:length(grid_sizes)
    disp([num2str(grid_sizes(i)), 'x', num2str(grid_sizes(i)), 'x', num2str(grid_sizes(i)), ...
          '    ', num2str(elapsed_times(i)), '           ', num2str(residual_norms(i)), ...
          '           ', num2str(loop_counts(i))]);
end
