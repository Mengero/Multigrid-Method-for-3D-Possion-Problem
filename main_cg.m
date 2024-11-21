%Nx = 500; Ny = 5; Nz = 5;  % anisotropic cell sizes
Nx = 32; Ny = 32; Nz = 32; % isotropic cell sizes

Ax = tridiag(Nx); Ix = speye(Nx);
Ay = tridiag(Ny); Iy = speye(Ny);
Az = tridiag(Nz); Iz = speye(Nz);

A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);

m = Nx * Ny * Nz;
f = rand(m, 1);
x0 = zeros(size(f));

[ret, iter] = cg(A_true,f,x0,1e-8,10000);

disp(norm(A_true * ret - f));