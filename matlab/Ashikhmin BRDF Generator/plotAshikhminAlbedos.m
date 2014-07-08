function [albedoes, thetas, phis] = ...
        plotAshikhminAlbedos(sigma, nGridTheta, nGridPhi, nTheta, nPhi)
    gGrid = Ashikhmin_gGrid(@(h) pGaussian(h, sigma), nGridTheta, nGridPhi);
    avgNH = AshikhminAverageNH(gGrid.densityFunction);
    
    albedoes = zeros(nTheta, nPhi);
    [phis, thetas] = meshgrid( ...
        linspace(0, pi/2, nPhi), linspace(0, pi/2, nTheta));
    wis = sph2vector(thetas, phis);
    for i=1:numel(thetas)
        fprintf('(%.2f, %.2f)\n', rad2deg(thetas(i)), rad2deg(phis(i)));
        albedoes(i) = AshikhminBSDFAlbedoUsingGrid(1.0, ...
            wis(:, i), gGrid, avgNH);
    end
    surf(rad2deg(thetas), rad2deg(phis), albedoes);
end