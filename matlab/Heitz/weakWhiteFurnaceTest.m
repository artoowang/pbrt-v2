% wo: 3x1 vector
function result = weakWhiteFurnaceTest(wo, D, G1)
    % Validity check
    if (wo(3) < 0)
        fprintf('wo should be in the upper hemisphere\n');
        return;
    end

    result = quad2d(@(theta, phi) integrand(theta, phi, wo, D, G1), ...
            0, pi/2, 0, 2*pi);
end

% wo: 3x1 vector
function val = integrand(theta, phi, wo, D, G1)
    wi = sph2vector(theta, phi); % 3xN
    wh = halfVector(wi, wo); % 3xN
    densities = D(wh); % 1xN
    G1_vals = G1(wo, wh, D); % 1xN
    
    % Reshape
    densities = reshape(densities, size(theta));
    G1_vals = reshape(G1_vals, size(theta));
    
    wg_dot_wo = wo(3);
    
    val = densities .* G1_vals ./ (4*wg_dot_wo) .* sin(theta);
end