function beta_list = beta_list_generation(nsmooth)
    beta_list = opt_coef(nsmooth);
end


function beta = opt_coef(k)
% Reference: Lottes, James. "Optimal polynomial smoothers for multigrid V‐cycles." 
%            Numerical Linear Algebra with Applications 30.6 (2023): e2518.
    x = cos([1:k]'*pi/(k+.5));
    w = (1-x)/(k+.5);
    
    W = zeros(k,k);
    W(:,1) = 1;
    if k >=2; W(:, 2) = 2*x + 1; end
    for i = 3:k
        W(:, i) = (2*x) .* W(:, i-1) - W(:, i-2);
    end

    r = opt_roots(k);
    lambda = (1-x)/2;
    p = prod(1-lambda'./r,1)';
    alpha = W' * (w .* p);
    beta = 1-cumsum((2*[0:k-1]'+1).*alpha, 1);

end

function r = opt_roots(k)
    function [w,f,g,ngp] = vars(r, x)
        p = prod(1- x' ./ r, 1)';
        w = x ./ (1- p.^2);
        f = sqrt(w) .* p;
        q = sum(1 ./ (x'- r), 1)';
        g = x .* (1./(2*w) + q);
        ngp = sum((p'.^2 + r./(x'-r))./(x'-r), 1)';
    end
    r = .5- .5*cos([1:k]'/(k+.5) * pi);
    x = .5- .5*cos((.5+[1:k-1]')/(k+.5) * pi);
    dr = r; drsize = 1;
    while drsize > 128*eps
        dx = x; dxsize = 1;
        while dxsize > 128*eps
            dxsize = norm(dx, inf);
            [~,~,g,ngp] = vars(r,x);
            dx = g ./ ngp; x = x + dx;
        end
        x1 = [x;1]; [w,f,~,~] = vars(r,x1);
        f0 = sqrt(.5/sum(1./r));
        J = f0^3./r'.^2 + w.*abs(f)./(r'.*(x1-r'));
        drsize = norm(dr, inf);
        dr =-J \ (f0- abs(f)); r = r + dr;
    end
end