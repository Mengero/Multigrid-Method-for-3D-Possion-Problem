function ret = multigridPreconditioner(A, Nx, Ny, Nz, f, mode, nsmooth)
    u = zeros(size(f));
    beta_list = beta_list_generation(nsmooth);
    [ret, ~] = vv_cycle(A, Nx, Ny, Nz, f, u, mode, nsmooth, beta_list);
end