% axes: 3xN matrix
% NDF: thetaRes x phiRes matrix
function NDF = gaussianNDFGenerator(thetaRes, phiRes, axes, sigmas)
    nAxes = size(axes, 2);
    NDF = zeros(1, thetaRes * phiRes);
    
    [pp, tt] = meshgrid(linspace(0, 2*pi, phiRes), linspace(0, pi, thetaRes));
    cellAxes = sph2vector(tt, pp);
    
    for i = 1:nAxes
        axis = axes(:, i);
        sigma = sigmas(i);
        NDF = NDF + pGaussian(cellAxes, sigma, axis);
    end
    
    NDF = reshape(NDF, [thetaRes phiRes]);
end