function vals = testAshikhminP(sigma, nTheta, nPhi)
    % Create a grid for P
    thetas = linspace(0, pi/2, nTheta);
    phis = linspace(0, 2*pi, nPhi);
    vals = zeros(nTheta, nPhi);
    for y=1:length(thetas)
        theta = thetas(y);
        for x=1:length(phis)
            phi = phis(x);
            k = sph2vector(theta, phi);
            vals(y, x) = AshikhminP(k, @(h) pGaussian(h, sigma));
        end
    end
    [pp, tt] = meshgrid(phis, thetas);
    surf(tt, pp, vals);
end