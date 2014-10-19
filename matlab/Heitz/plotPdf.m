function img = plotPdf(wo, Pdf, thetaRes, phiRes)
    %thetas = linspace(0, pi, thetaRes+1);
    %thetas = 0.5 * (thetas(1:end-1) + thetas(2:end));
    %phis = linspace(0, 2*pi, phiRes+1);
    %phis = 0.5 * (phis(1:end-1) + phis(2:end));
    thetas = linspace(0, pi, thetaRes);
    phis = linspace(0, 2*pi, phiRes);
    [pp, tt] = meshgrid(phis, thetas);
    wis = sph2vector(tt, pp);
    pdfs = Pdf(wo, wis);
    
    % Normalize
    %img = pdfs / sum(pdfs(:));
    %img = reshape(img, size(tt));
    %imagesc(img);
    
    x = reshape(wis(1,:), size(tt));
    y = reshape(wis(2,:), size(tt));
    z = reshape(wis(3,:), size(tt));
    img = reshape(pdfs, size(tt));
    surf(x, y, z, img, 'EdgeColor', 'none'), ...
        axis square, ...
        xlabel('X'), ylabel('Y'), zlabel('Z');
end