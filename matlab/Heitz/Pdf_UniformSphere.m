% wi: 3xN vectors
% pdf: 1xN vector
function pdf = Pdf_UniformSphere(wi)
    pdf = repmat(1/(4*pi), [1 size(wi, 2)]);
end