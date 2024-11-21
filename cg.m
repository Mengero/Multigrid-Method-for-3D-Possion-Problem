function [x, k] = cg(A, b, x0, tol, max_iter)
    x = x0;
    r = b - A * x;
    p = r;
    k = 0;

    while norm(r) > tol && k < max_iter
        alpha = (r' * r) / (p' * (A * p));
        x = x + alpha * p;
        r_new = r - alpha * (A * p);
        if norm(r_new) < tol
            break;
        end
        beta = (r_new' * r_new) / (r' * r);
        p = r_new + beta * p;
        r = r_new;
        k = k + 1;
    end
    if k == max_iter
        disp('Maximum iterations reached without convergence.');
    else
        disp(['Converged in ', num2str(k), ' iterations.']);
    end
end