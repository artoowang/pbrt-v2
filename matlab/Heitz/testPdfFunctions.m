function vals = testPdfFunctions(name, thetaRes, params)
    thetas = linspace(0, pi/2, thetaRes+1);
    thetas = 0.5 * (thetas(1:end-1) + thetas(2:end));
    
    % GGX
    if (strcmp(name, 'GGX'))
        alpha = params(1);
        vals = runPdfTest(thetas, ...
            @(wo, wi) Pdf_GGX(wo, wi, alpha), ...
            @(wh) D_GGX(wh, alpha), ...
            sprintf('GGX (alpha = %f)', alpha));
    end
end

function vals = runPdfTest(thetas, Pdf, D, name)
    thetaRes = length(thetas);
    vals = zeros(1, thetaRes);
    for t = 1:thetaRes
        fprintf('%d/%d\n', t, thetaRes);
        wo = sph2vector(thetas(t), 0);
        vals(t) = quad2d(...
            @(theta, phi) integrand(theta, phi, wo, Pdf), ...
            0, pi, 0, 2*pi);
        vals(t) = vals(t) + quad2d(...
            @(theta, phi) integrand2(theta, phi, wo, D), ...
            0, pi, 0, 2*pi);
    end
    plot(thetas, vals), ...
        xlabel('theta'), ...
        title(sprintf('Pdf integration result over wo (%s)', name)), ...
        ylim([0.9 1.1]);
end

% This integrand computes the differential probability mass of given wo and
% wi
function val = integrand(theta, phi, wo, Pdf)
    wi = sph2vector(theta, phi); % 3xN
    pdfs = Pdf(wo, wi); % 1xN
    
    % Reshape
    pdfs = reshape(pdfs, size(theta));
    
    val = pdfs .* sin(theta);
end

% This integrand computes the differential probability mass of given wh.
% Note: we only computes the wh facing away from wo
% Note 2: this indicator function creates a discontinuity on the function,
% which seems to cause problem for the quad2d near grazing angle
function val = integrand2(theta, phi, wo, D)
    wh = sph2vector(theta, phi); % 3xN
    densities = D(wh); % 1xN
    wo_dot_wh = wo' * wh;
    densities(wo_dot_wh > 0) = 0;
    
    % Reshape
    densities = reshape(densities, size(theta));
    
    val = densities .* cos(theta) .* sin(theta);
end