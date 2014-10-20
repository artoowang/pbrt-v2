function pdfs = plotPdf(wo, Pdf, thetaRes, phiRes)
    thetas = linspace(0, pi, thetaRes+1);
    thetas = 0.5 * (thetas(1:end-1) + thetas(2:end));
    phis = linspace(0, 2*pi, phiRes+1);
    phis = 0.5 * (phis(1:end-1) + phis(2:end));
    %thetas = linspace(0, pi, thetaRes);
    %phis = linspace(0, 2*pi, phiRes);
    [pp, tt] = meshgrid(phis, thetas);
    wis = sph2vector(tt, pp);
    pdfs = Pdf(wo, wis);
    pdfs = reshape(pdfs, size(tt));
    
    %imagesc(pdfs);
    
    %
    % Pad additional row and column to pdfs
    % Note: I think these padded values are not used by surf, but they are
    %       necessary to match the dimension of tt and pp
    %pdfs(:, end+1) = pdfs(:, end);
    %pdfs(end+1, :) = pdfs(end, :);
    pdfs(:, end+1) = pdfs(:, end);
    pdfs(end+1, :) = pdfs(end, :);
    
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
    surf(x, y, z, pdfs, 'EdgeColor', 'none'), ...
        axis square, ...
        xlabel('X'), ylabel('Y'), zlabel('Z');
    %}
end