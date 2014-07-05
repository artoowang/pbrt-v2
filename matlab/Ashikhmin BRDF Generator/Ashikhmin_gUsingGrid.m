% k: 3xN matrix
% gVal: 1xN vector
function gVal = Ashikhmin_gUsingGrid(k, gGrid)
    [thetas, phis] = vector2sph(k); % thetas, phis: 1xN vectors
    gVal = interp2(gGrid.phis, gGrid.thetas, gGrid.vals, phis, thetas);
end