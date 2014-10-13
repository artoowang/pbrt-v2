function img = plotG1 (wn, D, thetaRes, phiRes)
    thetas = linspace(0, pi/2, thetaRes+1);
    thetas = 0.5 * (thetas(1:end-1) + thetas(2:end));
    phis = linspace(0, 2*pi, phiRes+1);
    phis = 0.5 * (phis(1:end-1) + phis(2:end));
    %[pp, tt] = meshgrid(phis, thetas);
    img = zeros(thetaRes, phiRes);
    
    for t = 1:thetaRes
        for p = 1:phiRes
            wo = sph2vector(thetas(t), phis(p));
            img(t, p) = G1(wo, wn, D);
        end
    end
    
    %imagesc(img);
    %imshow(img);
end