function vals = testAshikhmin_g(sigma, nTheta, nPhi)
    % Create a image for g
    thetas = linspace(0, pi, nTheta);
    phis = linspace(0, 2*pi, nPhi);
    vals = zeros(nTheta, nPhi);
    for y=1:length(thetas)
        theta = thetas(y);
        for x=1:length(phis)
            phi = phis(x);
            k = sph2vector(theta, phi);
            vals(y, x) = Ashikhmin_g(k, @(h) pGaussian(h, sigma));
        end
    end
    [pp, tt] = meshgrid(phis, thetas);
    surf(tt, pp, vals);
end