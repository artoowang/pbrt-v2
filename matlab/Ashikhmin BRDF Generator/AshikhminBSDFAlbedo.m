function albedo = AshikhminBSDFAlbedo(sigma, wr)
    % Create a grid for P
    albedo = quad2d(@(theta, phi) integrand(theta, phi, sigma, wr), ...
        0, pi/2, 0, 2*pi);
end

function val = integrand(theta, phi, sigma, wr)
    val = zeros(size(theta));
    for i=1:numel(theta)
        wi = sph2vector(theta(i), phi(i));
        bsdfVal = AshikhminBSDF(wi, wr, @(h) pGaussian(h, sigma), 1);
        val(i) = bsdfVal * cos(theta(i));
    end
end