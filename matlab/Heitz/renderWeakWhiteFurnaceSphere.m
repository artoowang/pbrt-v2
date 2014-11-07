% f: @(wi, wo)
% Sample: @(wo, N)
function vals = renderWeakWhiteFurnaceSphere(f, Sample, resolution, N)
    thetas = linspace(0, pi/2, resolution+1);
    thetas = (thetas(1:end-1) + thetas(2:end)) / 2;
    wos = sph2vector(thetas, 0);
    
    vals = zeros(1, resolution);
    for i = 1:resolution
        fprintf('%d/%d ...\n', i, resolution);
        
        wo = wos(:, i);
        [wis, pdfs] = Sample(wo, N);
        wis = wis(:, pdfs > 0);
        pdfs = pdfs(pdfs > 0);
        if (isempty(pdfs))
            error('Monte Carlo doesn''t find any valid sample\n');
        end
        f_vals = f(wis, wo);  % 1xN
        costhetai = abs(wis(3, :));  % 1xN
        vals(i) = sum(f_vals .* costhetai ./ pdfs) / N;
    end
    
    drawCircle(0.5 * vals);
end