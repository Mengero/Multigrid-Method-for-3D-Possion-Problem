function [x, iter, relres] = myPCG(A, b, tol, max_iter, precond)
    x = zeros(size(b)); r = b; z = precond(r);
    p = z; rz_old = r' * z;
    
    for iter = 1:max_iter
        Ap = A * p;
        alpha = rz_old / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;

        relres = norm(r) / norm(b);
        if relres < tol
            disp(['Converged at iteration ', num2str(iter)]);
            return;
        end
        
        z = precond(r);
        rz_new = r' * z;
        beta = rz_new / rz_old;
        p = z + beta * p;
        rz_old = rz_new;
    end

    disp('Maximum iterations reached');
end
