function [ret, residue] = v_cycle(Ax, Ay, Az, f, u, mode, beta_list)
    % mode = 0 damped-Jacobi
    % mode = 1 Chebyshev Smoothing

    [Nx, ~] = size(Ax); [Ny, ~] = size(Ay); [Nz, ~] = size(Az);  
    Ix = speye(Nx); Iy = speye(Ny); Iz = speye(Nz);
    A = kron(kron(Az, Iy), Ix) + kron(kron(Iz, Ay), Ix) + kron(kron(Iz, Iy), Ax);
    
    if Nx == 1
        ret = A \ f;
        return; 
    end

    f_0 = A * u;

    if mode == 0
        r = f-f_0;
        omega = 2/3;
        d = diag(A);
        d_inv = 1 ./ d;
    elseif mode == 1
        d = diag(A); di=1./d;
        lmin = 0.4 * 3; lmax = 1.9 * 3;
        theta = .5*(lmax+lmin);
        delta = .5*(lmax-lmin);
        sigma = theta/delta;
        rho1  = 1./sigma;
        r=di.*(f-f_0); d=r/theta;
    elseif mode == 2
        r = f-f_0;
        lmax = 3*1.9; di = 1./diag(A);
        d = 4/3/lmax*di.*r;
        beta = 1; % 4th-kind without opt.
    end

    % interpolation matrix
    % Nx, Ny and Nz are the number of points inside the grid (not count the
    % boundary
        
    xc = 1 * [1:(Nx/2)]' / (Nx/2 + 1); x = 1 * [1:Nx]' / (Nx + 1); Jx = lin_interp_mat(xc,x);
    yc = 1 * [1:(Ny/2)]' / (Ny/2 + 1); y = 1 * [1:Ny]' / (Ny + 1); Jy = lin_interp_mat(yc,y);
    zc = 1 * [1:(Nz/2)]' / (Nz/2 + 1); z = 1 * [1:Nz]' / (Nz + 1); Jz = lin_interp_mat(zc,z);

    J = kron(Jz, kron(Jy, Jx));
    nsmooth=5; 
    
    
    % Pre-smoothing
    if mode == 0
        for nu=1:nsmooth
            u=u+omega*(d_inv .* r);
            r=f-A*u;
        end
    elseif mode == 1
        for k=1:800
            u = u + d;
            r = r - di.*(A*d);  rho0=rho1; rho1 = 1/(2*sigma-rho0);
            d = (rho1*rho0)*d + (2*rho1/delta)*r;
        end
        u = u+d;
    elseif mode == 2
        for k=1:2000
            % beta = beta_list(k);
            u = u + beta*d;
            r = r - A*d;
            display(norm(u))
            d = (2*k-1)/(2*k+3)*d + (8*k+4)/(2*k+3)/lmax*di.*r;
        end
        u = u + beta*d;
    end
    % Restriction
    rhs = J * r; 
    eps = zeros(size(rhs));
    
    if Nx == 1
        % assume Ny == Nz here
        Ac = J * A * J';
        eps = Ac \ temp;
    else
        AAx = tridiag(Nx/2); AAy = tridiag(Ny/2); AAz = tridiag(Nz/2);
        eeps = zeros(size(eps));
        eps = v_cycle(AAx, AAy, AAz, rhs, eeps, mode, beta_list);
    end
    ec = J' * eps;
    u = u + ec;
    
    if mode == 0
        r = f - A * u;
    elif 
        r = di.*(f - A * u);
    end
    
    % Post-Smoothing
    if mode == 0
        for nu=1:nsmooth
            u=u+omega*(d_inv .* r);
            r=f-A*u;
        end
    elseif mode == 1
        for k=1:800
            u = u + d;
            r = r - di.*(A*d);  rho0=rho1; rho1 = 1/(2*sigma-rho0);
            d = (rho1*rho0)*d + (2*rho1/delta)*r;
        end
        u = u+d;
    elseif mode == 2
        for k=1:2000
            % beta = beta_list(k);
            u = u + beta*d;
            r = r - A*d;
            d = (2*k-1)/(2*k+3)*d + (8*k+4)/(2*k+3)/lmax*di.*r;
        end
        u = u + beta*d;
    end
    
    ret = u;
    residue = r;
end