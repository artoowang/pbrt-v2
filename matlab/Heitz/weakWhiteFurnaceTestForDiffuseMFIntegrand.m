% wo: 3x1 vector
% wis: 3xN vectors
function val = weakWhiteFurnaceTestForDiffuseMFIntegrand(wis, wo, D, G1)
    N = size(wis, 2);
    val = zeros(1, N);
    
    wo_dot_wg_abs = abs(wo(3));
    if (wo_dot_wg_abs > 0)
        for i = 1:N
            val(i) = quad2d(@(theta, phi) integrandOverWm(theta, phi, wis(:, i), wo, D, G1), ...
                0, pi, 0, 2*pi);
        end
        val = val / (pi * wo_dot_wg_abs);
    end
end

function val = integrandOverWm(theta, phi, wi, wo, D, G1)
    wm = sph2vector(theta, phi);
    
    wo_dot_wm_clamped = wo' * wm;  % 1xN
    wo_dot_wm_clamped(wo_dot_wm_clamped < 0) = 0;
    wi_dot_wm_clamped = wi' * wm;  % 1xN
    wi_dot_wm_clamped(wi_dot_wm_clamped < 0) = 0;
    
    G1_vals = G1(wo, wm, D);
    D_vals = D(wm);
    
    val = wo_dot_wm_clamped .* wi_dot_wm_clamped .* G1_vals .* D_vals .* sin(theta(:)');
    val = reshape(val, size(theta));
end