function testDFunctions(name, params)
    thetas = linspace(0, pi, 100);
    
    % Beckmann
    if (strcmp(name, 'Beckmann'))
        alpha = params(1);
        ws = sph2vector(thetas, zeros(size(thetas)));
        plot(thetas, D_Beckmann(ws, alpha)), xlabel('theta'), title(sprintf('Beckmann (alpha = %f)', alpha));
        val = quad2d(@(theta, phi) integrand(theta, phi, @(w) D_Beckmann(w, alpha)), ...
            0, pi, 0, 2*pi);
        fprintf('Beckmann (alpha = %f) integrates to %f\n', alpha, val);
    end
    
    % GGX
    if (strcmp(name, 'GGX'))
        alpha = params(1);
        ws = sph2vector(thetas, zeros(size(thetas)));
        plot(thetas, D_GGX(ws, alpha)), xlabel('theta'), title(sprintf('GGX (alpha = %f)', alpha));
        val = quad2d(@(theta, phi) integrand(theta, phi, @(w) D_GGX(w, alpha)), ...
            0, pi, 0, 2*pi);
        fprintf('GGX (alpha = %f) integrates to %f\n', alpha, val);
    end

    % Blinn
    if (strcmp(name, 'Blinn'))
        exponent = params(1);
        ws = sph2vector(thetas, zeros(size(thetas)));
        plot(thetas, D_Blinn(ws, exponent)), xlabel('theta'), title(sprintf('Blinn (exponent = %f)', exponent));
        val = quad2d(@(theta, phi) integrand(theta, phi, @(w) D_Blinn(w, exponent)), ...
            0, pi, 0, 2*pi);
        fprintf('Blinn (exponent = %f) integrates to %f\n', exponent, val);
    end
end

function val = integrand(theta, phi, D)
    w = sph2vector(theta, phi); % 3xN
    densities = D(w); % 1xN
    
    % Reshape
    densities = reshape(densities, size(theta));
    
    val = densities .* cos(theta) .* sin(theta);
end