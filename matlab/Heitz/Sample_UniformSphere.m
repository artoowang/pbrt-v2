% wo: 3x1 vector
% wi: 3xN vectors
% pdf: 1xN vector
function [wi, pdf] = Sample_UniformSphere(N)
    u1 = rand(1, N);
    u2 = rand(1, N);
    z = 1 - 2 * u1;
    r = sqrt(max(0, 1 - z .* z));
    phi = 2 * pi * u2;
    wi = [r .* cos(phi); r .* sin(phi); z];
    pdf = repmat(1/(4*pi), [1 N]);
end