%Nx = 500; Ny = 5; Nz = 5;  % anisotropic cell sizes
Nx = 32; Ny = 32; Nz = 32; % isotropic cell sizes

Ax = tridiag(Nx); Ix = speye(Nx);
Ay = tridiag(Ny); Iy = speye(Ny);
Az = tridiag(Nz); Iz = speye(Nz);

m = Nx * Ny * Nz;
f = rand(m, 1);

[u, time] = FDM(Ax, Ay, Az,f);
disp(time);
[u, time] = FDM(Ax, Ay, Az,f);
disp(time); % use this after the warming

% verify
A_true = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);

disp(norm(u - A_true \ f));
