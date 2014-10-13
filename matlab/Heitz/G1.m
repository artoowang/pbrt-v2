% wo, wn is assumed to be in the tangent frame (where surface normal wg is +Z)
% wo: 3x1
% wn: 3xN
% Note in general G1 does not depend on wn (i.e., the normal/masking
% independence assumption); wn is only used to cull backfacing facets
function results = G1(wo, wn, D)
    denominator = quad2d(@(theta, phi) integrand(theta, phi, wo, D), ...
        0, pi/2, 0, 2*pi);
    
    if (denominator < 1e-6)
        fprintf('denominator is too small: %e\n', denominator);
    end
    
    cos_theta_o = wo(3);
    wo_dot_wn = wo' * wn;
    results = zeros(1, size(wn, 2));
    results(wo_dot_wn > 0) = cos_theta_o ./ denominator;
    
    if (any(isnan(results)))
        fprintf('Error: G1 has NaN\n');
    end
end

% wo: 3x1 vector
function val = integrand(theta, phi, wo, D)
    wn = sph2vector(theta, phi); % 3xN
    wo_dot_wn_clamped = wo' * wn; % 1xN
    wo_dot_wn_clamped(wo_dot_wn_clamped < 0) = 0;
    densities = D(wn); % 1xN
    
    % Reshape
    wo_dot_wn_clamped = reshape(wo_dot_wn_clamped, size(theta));
    densities = reshape(densities, size(theta));
    
    val = wo_dot_wn_clamped .* densities .* sin(theta);
end