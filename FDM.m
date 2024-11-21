function [ret, time] = FDM(Ax, Ay, Az, f)
    tic;
    [Nx, ~] = size(Ax);
    [Ny, ~] = size(Ay);
    [Nz, ~] = size(Az);
    
    m = Nx * Ny * Nz;
    
    Ax_full=full(Ax); [Sx,Lamx]=eig(Ax_full); Lamx = diag(Lamx); Sxi=inv(Sx);
    Ay_full=full(Ay); [Sy,Lamy]=eig(Ay_full); Lamy = diag(Lamy); Syi=inv(Sy);
    Az_full=full(Az); [Sz,Lamz]=eig(Az_full); Lamz = diag(Lamz); Szi=inv(Sz);
    
    test_A = kron(Szi, kron(Syi, Sxi));
    test_f = test_A * f;

    f = Kronecker_3D(Sxi, Syi, Szi, Nx, Ny, Nz, f);

    D = sparse(m, m);
    for j = 1:Nz
        for k = 1:Ny
            for p = 1:Nx
                index = (j-1)*Nx*Ny+(k-1)*Nx+p;
                D(index,index) = 1/(Lamx(p) + Lamy(k) + Lamz(j));
            end
        end
    end
    
    f = D * f;

    ret = Kronecker_3D(Sx, Sy, Sz, Nx, Ny, Nz, f);
    time = toc;
end