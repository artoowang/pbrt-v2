% wo: 3x1 vector
% wis: 3xN vectors
% D: @(w)
% G1: @(wo, wn, D)
% F: @(costheta)
% vals: 1xN vector
function val = f_G1_diffuse(wis, wo, D, G1, F)
    N = size(wis, 2);
    val = zeros(1, N);
    
    wo_dot_wg_abs = abs(wo(3));
    if (wo_dot_wg_abs > 0)
        wis_dot_wg_abs = abs(wis(3, :));
        one_over_wis_dot_wg_abs = 1 ./ wis_dot_wg_abs;
        one_over_wis_dot_wg_abs(wis_dot_wg_abs == 0) = 0;
        
        for i = 1:N
            % Method 1: quad2d
            %
            val(i) = quad2d(@(theta, phi) integrandOverWmThetaPhi(...
                theta, phi, wis(:, i), wo, D, G1, F), ...
                0, pi, 0, 2*pi);
            %}
            % Method 2: Monte Carlo
            %   Note this can use either uniform sphere sampling (very
            %   inefficient) or Sample_D_GGX, but the later actually needs
            %   the alpha value (roughness), which is not passed into this
            %   function so is therefore hard coded
            %{
            val(i) = monteCarloIntegrate(...
                @(w) integrandOverWm(w, wis(:, i), wo, D, G1, F), ...
                ... %@Sample_UniformSphere,
                @(N) Sample_D_GGX(N, 0.1), ...
                1000);
            %}
        end
        
        val = val .* one_over_wis_dot_wg_abs / (pi * wo_dot_wg_abs);
    end
end

% wo, wi: 3x1 vector
function val = integrandOverWm(wm, wi, wo, D, G1, F)
    wo_dot_wm = wo' * wm;  % 1xN
    wo_dot_wm_clamped = wo_dot_wm;
    wo_dot_wm_clamped(wo_dot_wm < 0) = 0;
    wi_dot_wm = wi' * wm;  % 1xN
    wi_dot_wm_clamped = wi_dot_wm;
    wi_dot_wm_clamped(wi_dot_wm < 0) = 0;
    
    G1_vals = G1(wo, wm, D);
    D_vals = D(wm);
    F_vals = F(wi_dot_wm);
    
    val = wo_dot_wm_clamped .* wi_dot_wm_clamped .* G1_vals .* D_vals .* F_vals;
end

function val = integrandOverWmThetaPhi(theta, phi, wi, wo, D, G1, F)
    wm = sph2vector(theta, phi);
    val = integrandOverWm(wm, wi, wo, D, G1, F) .* sin(theta(:)');
    val = reshape(val, size(theta));
end