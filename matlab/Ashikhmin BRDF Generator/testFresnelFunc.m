function testFresnelFunc(sigma)
    f0 = 0.04;
    theta_i = linspace(0, pi/2, 20);
    wi_dot_N = cos(theta_i);
    plot(theta_i, Fresnel(wi_dot_N, f0)), xlabel('theta\_i');
    hold on;
    
    plot(theta_i, 1-Fresnel(wi_dot_N, f0), 'r-');
    
    % TODO: test
    %{
    intPGaussian = quad2d( ...
        @(theta, phi) integralPGaussian(theta, phi, sigma), ...
        0, pi/2, 0, 2*pi);
    fprintf('pGaussian(%.2f) integrates (over hemisphere) to %f\n', sigma, intPGaussian);
    %}
    
    %
    Ft = zeros(size(theta_i));
    wis = sph2vector(theta_i, zeros(size(theta_i)));
    for i = 1:length(theta_i)
    	Ft(i) = quad2d( ...
            @(theta, phi) integralFtTimesPGaussian(theta, phi, sigma, wis(:,i), f0), ...
            0, pi/2, 0, 2*pi);
    end
    plot(theta_i, Ft, 'g-');
    %}
    
    legend('Fr(dot(wi,N))', '1-Fr(dot(wi,N))', sprintf('Ft(sigma=%f)', sigma));
    hold off;
end

function val = integralPGaussian(theta, phi, sigma)
    v = sph2vector(theta, phi);
    densities = pGaussian(v, sigma);
    val = reshape(densities(:), size(theta)) .* sin(theta);
end

% wi: 3x1
function val = integralFtTimesPGaussian(theta, phi, sigma, wi, f0)
    wh = sph2vector(theta, phi);    % 3xN
    densities = pGaussian(wh, sigma);  % 1xN
    wi_dot_wh = wi' * wh; % 1xN
    Fr = Fresnel(wi_dot_wh, f0); % 1xN
    tmp = (1-Fr).*densities;
    % TODO: test
    tmp(isnan(tmp)) = 0;
    if (sum(isnan(tmp(:))) > 0)
        fprintf('Error: NaN\n');
    end
    val = reshape(tmp, size(theta)) .* sin(theta);
end