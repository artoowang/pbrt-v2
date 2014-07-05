function testPGaussian(sigma)
    % Computes integration
    integration = quad2d( ...
        @(theta, phi) integral(theta, phi, sigma), ...
        0, pi/2, 0, 2*pi);
    fprintf('integrates to %f over hemisphere\n', integration);
    
    % Show value on theta/phi image
    theta = linspace(0, pi, 512);
    phi = linspace(0, 2*pi, 1024);
    [tt, pp] = meshgrid(theta, phi);
    [x, y, z] = sph2cart(pp(:), pi/2-tt(:), 1);
    h = [x'; y'; z'];
    densities = pGaussian(h, sigma);
    densities = reshape(densities, size(tt));
    imagesc(densities);
end

function val = integral(theta, phi, sigma)
    [x, y, z] = sph2cart(phi(:), pi/2-theta(:), 1);
    h = [x'; y'; z'];
    densities = pGaussian(h, sigma);
    val = reshape(densities(:), size(theta)) .* sin(theta);
end