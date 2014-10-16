% wo: 3x1 vector
function result = weakWhiteFurnaceTest(wo, D, G1, thetaRes, phiRes)
    % Validity check
    if (wo(3) < 0)
        fprintf('wo should be in the upper hemisphere\n');
        return;
    end

    % Method 1: quad2d
    %{
    [result, errbnd] = quad2d(@(theta, phi) integrand(theta, phi, wo, D, G1), ...
            ...%1, 2, 3, 3.5);
            0, pi, 0, 2*pi);
    %errbnd
    %}
        
    % Method 2: fixed grid
    %
    if (~exist('thetaRes', 'var'))
        thetaRes = 200;
    end
    if (~exist('phiRes', 'var'))
        phiRes = 200;
    end
    thetas = linspace(0, pi, thetaRes+1);
    thetas = (thetas(1:end-1) + thetas(2:end)) / 2;
    %dtheta = 0.05;
    %thetas = 0:dtheta:pi;
    phis = linspace(0, 2*pi, phiRes+1);
    phis = (phis(1:end-1) + phis(2:end)) / 2;
    %dphi = 0.05;
    %phis = 0:dphi:2*pi;
    [pp, tt] = meshgrid(phis, thetas);
    wis = sph2vector(tt, pp);
    whs = halfVector(wis, wo);
    densities = D(whs);
    G1_vals = G1(wo, whs, D);
    wg_dot_wo = wo(3);
    vals = densities .* G1_vals ./ (4*wg_dot_wo) .* sin(tt(:)');
    %surf(pp, tt, reshape(vals, size(tt)));
    result = sum(vals(:)) * pi / thetaRes * 2*pi / phiRes;
    %result = sum(vals(:)) * dphi * dtheta;
    %}
end

% wo: 3x1 vector
function val = integrand(theta, phi, wo, D, G1)
    wg_dot_wo = wo(3);
    if (wg_dot_wo < 0)
        val = zeros(size(theta));
        return;
    end
    
    wi = sph2vector(theta, phi); % 3xN
    wh = halfVector(wi, wo); % 3xN
    densities = D(wh); % 1xN
    G1_vals = G1(wo, wh, D); % 1xN
    
    % Reshape
    densities = reshape(densities, size(theta));
    G1_vals = reshape(G1_vals, size(theta));
    
    val = densities .* G1_vals ./ (4*wg_dot_wo) .* sin(theta);
end