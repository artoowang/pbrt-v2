function albedo = AshikhminBSDFAlbedo(sigma, f0, wr, nTheta, nPhi)
    % Create one-time data structures
    gGrid = Ashikhmin_gGrid(@(h) pGaussian(h, sigma), nTheta, nPhi);
    avgNH = AshikhminAverageNH(gGrid.densityFunction);
    albedo = AshikhminBSDFAlbedoUsingGrid(f0, wr, gGrid, avgNH);
end