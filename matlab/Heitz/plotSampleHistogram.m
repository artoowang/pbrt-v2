function histmat = plotSampleHistogram(wis, thetaRes, phiRes)
    
    N = size(wis, 2);
    [theta_is, phi_is] = vector2sph(wis);
    histmat = hist2(phi_is, theta_is, [0 2*pi], [0 pi], phiRes, thetaRes);
    
    % Adjustment
    thetas = linspace(0, pi, thetaRes+1);
    thetas = (thetas(1:end-1) + thetas(2:end)) / 2;
    phis = linspace(0, 2*pi, phiRes+1);
    phis = (phis(1:end-1) + phis(2:end)) / 2;
    dphi = (2*pi) / phiRes;
    dtheta = pi / thetaRes;
    dw = sin(thetas) * dphi * dtheta;
    probMass = histmat / N;
    probDensities = probMass ./ repmat(dw(:), [1 phiRes]);
    histmat = probDensities;
    
    imagesc(histmat);
    
    %{
    % Pad additional row and column to histmat
    % Note: I think these padded values are not used by surf, but they are
    %       necessary to match the dimension of tt and pp
    histmat(:, end+1) = histmat(:, end);
    histmat(end+1, :) = histmat(end, :);
    
    % Compute xyz points for surf()
    % Note these points represent the endpoints of the cell, not the center
    % we use above to adjust the histogram
    thetas = linspace(0, pi, thetaRes+1);
    phis = linspace(0, 2*pi, phiRes+1);
    [pp, tt] = meshgrid(phis, thetas);
    ws = sph2vector(tt, pp);
    x = reshape(ws(1,:), size(tt));
    y = reshape(ws(2,:), size(tt));
    z = reshape(ws(3,:), size(tt));
    surf(x, y, z, histmat, 'EdgeColor', 'none'), ...
        axis square, ...
        xlabel('X'), ylabel('Y'), zlabel('Z');
    %}
end