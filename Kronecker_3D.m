function ret = Kronecker_3D(Sx ,Sy, Sz, mx ,my, mz, u)
    v = zeros(size(u));
    u_reshaped = reshape(u, [mx, my, mz]);
    for k = 1:mz
        for j = 1:my
            u_block = u_reshaped(:, j, k);
            v_block = Sx * u_block;
            v((k-1)*mx*my + (j-1)*mx + 1 : (k-1)*mx*my + j*mx) = v_block;
        end
    end
    
    v_reshaped = reshape(v, [mx, my, mz]);
    w_reshaped = zeros(size(v_reshaped));

    for k = 1:mz
        for i = 1:mx
            v_block = v_reshaped(i, :, k)';
            w_slice = Sy * v_block;
            w_reshaped(i, :, k) = w_slice';
        end
    end
    z_reshaped = zeros(size(w_reshaped));
    for j = 1:my
        for i = 1:mx
            w_slice = squeeze(w_reshaped(i, j, :));
            z_slice = Sz * w_slice;
            z_reshaped(i, j, :) = z_slice; 
        end
    end
    
    ret = z_reshaped(:);
    % for i = 1:mz
    %     for j = 1:my
    %         for k = 1:mx
    %             for p = 1:mz
    %                 v((i-1)*mx*my+(j-1)*mx+k) = v((i-1)*mx*my+(j-1)*mx+k) + Sx(i,p) * u((p-1)*mx*my+(j-1)*mx+k);
    %             end
    %         end
    %     end
    % end
    % for i = 1:mz
    %     for j = 1:my
    %         for k = 1:mx
    %             for p = 1:my
    %                 w((i-1)*mx*my+(j-1)*mx+k) = w((i-1)*mx*my+(j-1)*mx+k) + Sy(j,p) * v((i-1)*mx*my+(p-1)*mx+k);
    %             end
    %         end
    %     end
    % end
    % for i = 1:mz
    %     for j = 1:my
    %         for k = 1:mx
    %             for p = 1:mx
    %                 ret((i-1)*mx*my+(j-1)*mx+k) = ret((i-1)*mx*my+(j-1)*mx+k) + Sz(k,p) * w((i-1)*mx*my+(j-1)*mx+p);
    %             end
    %         end
    %     end
    % end
end