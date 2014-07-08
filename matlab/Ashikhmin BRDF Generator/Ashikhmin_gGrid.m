function gGrid = Ashikhmin_gGrid(densityFunction, nTheta, nPhi)
    % Create a grid for g
    % Note: it seems there is no need to test theta over pi/2, since g() is
    %       only evaluated for either wi or wr, which are always in the
    %       upper hemisphere
    thetas = linspace(0, pi/2, nTheta);
    phis = linspace(0, 2*pi, nPhi);
    gGrid.vals = zeros(nTheta, nPhi);
    gGrid.densityFunction = densityFunction;
    
    [gGrid.phis, gGrid.thetas] = meshgrid(phis, thetas);
    N = numel(gGrid.thetas);
    k = sph2vector(gGrid.thetas(:), gGrid.phis(:));   % 3xN matrix
    for i=1:N
        gGrid.vals(i) = Ashikhmin_g(k(:, i), densityFunction);
    end
end