function histmat = plotSampleHistogram(wis, thetaRes, phiRes)
    thetas = linspace(0, pi, thetaRes);
    phis = linspace(0, 2*pi, phiRes);
    
    [theta_is, phi_is] = vector2sph(wis);
    histmat = hist2(theta_is, phi_is, thetas, phis);
    
    % Adjustment
    factors = 1 ./ sin(thetas);
    factors(thetas == 0) = 0;
    histmat = histmat .* repmat(factors(:), [1 phiRes]);
    
    % Normalize
    %img = pdfs / sum(pdfs(:));
    %img = reshape(img, size(tt));
    %imagesc(img);
    
    [pp, tt] = meshgrid(phis, thetas);
    ws = sph2vector(tt, pp);
    x = reshape(ws(1,:), size(tt));
    y = reshape(ws(2,:), size(tt));
    z = reshape(ws(3,:), size(tt));
    surf(x, y, z, histmat, 'EdgeColor', 'none'), ...
        axis square, ...
        xlabel('X'), ylabel('Y'), zlabel('Z');
end