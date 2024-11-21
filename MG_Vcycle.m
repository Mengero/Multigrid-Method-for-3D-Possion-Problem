function u = MG_Vcycle(level, A_list, P_list, b, u0, nsmooth, mode)
% Multigrid V-clcye
% INPUT:
    % level     : The current level (initial is 1)
    % A_list    : The array of coarsed A matrices on each level
    % P_list    : The array of interpolation (prolongation) operators on each level
    % b         : The right hand side
    % u0        : Initial guess
    % mode      : Smoothing mode
    % nsmooth   : Number of Iterations in smoothing

%OUTPUT:
    % u         : The new solution after a V-clcye

	% Load coefficient matrix
	A = cell2mat(A_list(level));

    % Pre-smoothing
    u = u0;
	if mode == 0
        omega = 2/3;
        d = diag(A);
        d_inv = 1 ./ d;
        for nu=1:nsmooth
            u=u+omega*(d_inv .* r);
            r=f-A*u;
        end
    else
        d = diag(A); di=1./d;
        lmax = 1.9;
        lmin = 0.4;
        theta = .5*(lmax+lmin);
        delta = .5*(lmax-lmin);
        sigma = theta/delta;
        rho1  = 1./sigma;
        for k=1:50
           u = u + d;
           r = r - di.*(A*d);  rho0=rho1; rho1 = 1/(2*sigma-rho0);
           d = (rho1*rho0)*d + (2*rho1/delta)*r;
        end
    end
	
	% Load restriction operator and construct interpolation operator
	P = cell2mat(P_list(level));
	R = P';
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	r   = b - A * u;
	r_C = R * r;
	
	% Solve coarse grid problem recursively
	u0  = zeros(coarse_n, 1);
	e_C = MG_Vcycle(level + 1, A_list, P_list, r_C, u0, nsmooth);
	
	% Transfer error to fine grid and correct
	u = u + P * e_C;
	
	% Post-smoothing
	if mode == 0
        for nu=1:nsmooth
            u=u+omega*(d_inv .* r);
            r=f-A*u;
        end
    else
        for k=1:50
           u = u + d;
           r = r - di.*(A*d);  rho0=rho1; rho1 = 1/(2*sigma-rho0);
           d = (rho1*rho0)*d + (2*rho1/delta)*r;
        end
    end
end