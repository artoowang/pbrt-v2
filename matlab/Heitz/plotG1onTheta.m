function vals = plotG1onTheta (wn, G1, D, thetaRes)
    thetas = linspace(0, pi/2, thetaRes+1);
    thetas = 0.5 * (thetas(1:end-1) + thetas(2:end));
    vals = zeros(1, thetaRes);
    
    for t = 1:thetaRes
        wo = sph2vector(thetas(t), 0);
        vals(t) = G1(wo, wn, D);
    end
    
    plot(thetas, vals), xlabel('theta'), title('G1');
end