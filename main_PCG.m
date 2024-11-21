Nx = 16; Ny = 64; Nz = 64;  % anisotropic cell sizes
% Nx = 32; Ny = 32; Nz = 32; % isotropic cell sizes

Ax = tridiag(Nx); Ix = speye(Nx);
Ay = tridiag(Ny); Iy = speye(Ny);
Az = tridiag(Nz); Iz = speye(Nz);

m = Nx * Ny * Nz;
f = rand(m, 1);
u = zeros(size(f));
A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);

mode = 0;

nsmooth = 5;

tic;
M = @(x) multigridPreconditioner(A_true, Nx, Ny, Nz, x, mode, nsmooth);
[x, iter, resl] = myPCG(A_true, f, 1e-8, 2000, M);
elapsed_time = toc;
disp(elapsed_time);

disp(['Solution converged in ', num2str(iter), ' iterations with relative residual ', num2str(resl)]);
disp(norm(A_true * x - f));

x0 = zeros(size(f));
[ret, iter] = cg(A_true,f,x0,1e-8,10000);
disp(iter);
disp(norm(A_true * ret - f));