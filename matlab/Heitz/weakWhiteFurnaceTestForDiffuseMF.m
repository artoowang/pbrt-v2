function result = weakWhiteFurnaceTestForDiffuseMF(wo, D, G1)
    [result, errbnd] = quad2d(@(theta, phi) integrand(theta, phi, wo, D, G1), ...
            0, pi, 0, 2*pi);
end

function val = integrand(theta, phi, wo, D, G1)
    wi = sph2vector(theta, phi);
    val = weakWhiteFurnaceTestForDiffuseMFIntegrand(wi, wo, D, G1);
    val = reshape(val, size(theta)) .* sin(theta);
end