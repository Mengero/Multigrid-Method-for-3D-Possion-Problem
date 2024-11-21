Nx = 16; Ny = 16; Nz = 16;  % anisotropic cell sizes
% Nx = 32; Ny = 32; Nz = 32; % isotropic cell sizes

Ax = tridiag(Nx); Ix = speye(Nx);
Ay = tridiag(Ny); Iy = speye(Ny);
Az = tridiag(Nz); Iz = speye(Nz);

m = Nx * Ny * Nz;
f = rand(m, 1);
u = zeros(size(f));
A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);

ret = zeros(size(f));
residue = f;
beta_list = beta_list_generation;

tic;
iteration = 0;
while norm(residue) > 1e-10
    iteration = iteration + 1;
    [ret, rr] = v_cycle(A_true, Nx, Ny, Nz, f, u, 1);
    u = ret;
    residue = rr;
end
elapsed_time = toc;
disp(elapsed_time);
disp(iteration)
disp(norm(A_true * ret - f));

% tic;
% while norm(A_true * ret - f) > 1e-8
%     [ret, time] = v_cycle(Ax, Ay, Az, f, u, 1);
%     u = ret;
%     % disp(time);
% end
% elapsed_time = toc;
% disp(elapsed_time);
% disp(norm(A_true * ret - f));