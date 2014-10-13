% wo, wn is assumed to be in the tangent frame (where surface normal wg is +Z)
% Currently we only support single 3x1 vectors wo and wo
function val = avgDotWoWn(wo, D)
    val = quad2d(@(theta, phi) integrand(theta, phi, wo, D), ...
        0, pi/2, 0, 2*pi);
end

% wo: 3x1 vector
function val = integrand(theta, phi, wo, D)
    w = sph2vector(theta, phi); % 3xN
    wo_dot_w = wo' * w; % 1xN
    densities = D(w); % 1xN
    
    % Reshape
    wo_dot_w = reshape(wo_dot_w, size(theta));
    densities = reshape(densities, size(theta));
    
    val = wo_dot_w .* densities .* sin(theta);
end