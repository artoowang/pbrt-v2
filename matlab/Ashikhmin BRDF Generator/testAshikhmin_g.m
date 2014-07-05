function vals = testAshikhmin_g(sigma, nTheta, nPhi)
    % Using Ashikhmin_g
    thetas = linspace(0, pi/2, nTheta);
    phis = linspace(0, 2*pi, nPhi);
    [pp, tt] = meshgrid(phis, thetas);
    N = numel(tt);
    k = sph2vector(tt(:), pp(:));   % 3xN matrix
    vals = zeros(nTheta, nPhi);
    for i=1:N
        vals(i) = Ashikhmin_g(k(:, i), @(h) pGaussian(h, sigma));
    end
    figure, surf(tt, pp, vals), title('Using Ashikhmin\_g');
    
    % Using Ashikhmin_gUsingGrid
    gGrid = Ashikhmin_gGrid(@(h) pGaussian(h, sigma), nTheta, nPhi);
    figure, surf(gGrid.thetas, gGrid.phis, gGrid.vals), title('gGrid');
    gVals = Ashikhmin_gUsingGrid(k, gGrid);
    figure, surf(tt, pp, reshape(gVals, size(tt))), title('Using Ashikhmin\_gUsingGrid');
end