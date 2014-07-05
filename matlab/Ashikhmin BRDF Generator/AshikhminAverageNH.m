function val = AshikhminAverageNH (densityFunction)
    val = quad2d(@(theta, phi) integrand(theta, phi, densityFunction), ...
        0, pi/2, 0, 2*pi);
end

function val = integrand(theta, phi, densityFunction)
    h = sph2vector(theta, phi);
    densities = densityFunction(h);
    densities = reshape(densities(:), size(theta));
    % N dot H is cos(theta) by definition
    val = cos(theta) .* densities .* sin(theta);
end