% axes: 3xN matrix
% NDF: thetaRes x phiRes matrix
function NDF = gaussianNDFGenerator(thetaRes, phiRes, axes, sigmas)
    nAxes = size(axes, 2);
    NDF = zeros(1, thetaRes * phiRes);
    
    phis = linspace(0, 2*pi, phiRes+1);
    phis = 0.5 * (phis(1:end-1) + phis(2:end));
    
    thetas = linspace(0, pi, thetaRes+1);
    thetas = 0.5 * (thetas(1:end-1) + thetas(2:end));
    
    [pp, tt] = meshgrid(phis, thetas);
    cellAxes = sph2vector(tt, pp);
    sinthetas = sin(tt(:)');
    
    for i = 1:nAxes
        axis = axes(:, i);
        sigma = sigmas(i);
        % TODO: do we need to multiply NDF by sin(theta)?
        %NDF = NDF + pGaussian(cellAxes, sigma, axis) .* sinthetas;
        NDF = NDF + pGaussian(cellAxes, sigma, axis);
    end
    
    NDF = reshape(NDF, [thetaRes phiRes]);
end