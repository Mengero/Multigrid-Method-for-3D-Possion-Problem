% Nx = 16; Ny = 16 * 4; Nz = 16 * 4; % anisotropic cell sizes
Nx = 32; Ny = 32; Nz = 32; % isotropic cell sizes
Ax = tridiag(Nx); Ix = speye(Nx);
Ay = tridiag(Ny); Iy = speye(Ny);
Az = tridiag(Nz); Iz = speye(Nz);
m = Nx * Ny * Nz;
f = rand(m, 1);
u = zeros(size(f));
A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);
ret = zeros(size(f));
residue = f;
ret = zeros(size(f));
tic;

mode = 1;
nsmooth = 50;

beta_list = beta_list_generation(nsmooth);

while norm(residue)/norm(f) > 1e-8
[u, rr] = vv_cycle(A_true, Nx, Ny, Nz, f, u, mode, nsmooth, beta_list);
residue = rr;
disp(norm(residue));
end
elapsed_time = toc;
disp(elapsed_time);
disp(norm(A_true * u - f));