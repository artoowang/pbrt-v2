function albedo = AshikhminBSDFAlbedoUsingGrid(f0, wr, gGrid, avgNH)
    albedo = quad2d(@(theta, phi) integrand(theta, phi, wr, gGrid, avgNH, f0), ...
        0, pi/2, 0, 2*pi);
end

function val = integrand(theta, phi, wr, gGrid, avgNH, f0)
    wi = sph2vector(theta, phi);
    bsdfVal = AshikhminBSDFUsingGrid(wi, wr, f0, gGrid, avgNH);
    val = reshape(bsdfVal, size(theta)) .* cos(theta) .* sin(theta);
end