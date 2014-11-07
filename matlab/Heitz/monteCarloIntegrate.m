% sphFunc: @(w)
% Sample: @(N)
function val = monteCarloIntegrate(sphFunc, Sample, N)
    [ws, pdfs] = Sample(N);
    ws = ws(:, pdfs > 0);
    pdfs = pdfs(pdfs > 0);
    if (isempty(pdfs))
        error('Monte Carlo doesn''t find any valid sample\n');
    end
    sphFuncVals = sphFunc(ws);  % 1xN
    val = sum(sphFuncVals ./ pdfs) / N;
end