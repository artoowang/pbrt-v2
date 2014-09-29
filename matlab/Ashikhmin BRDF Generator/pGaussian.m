% h_local is in tangent space
% (i.e., shading normal N is (0,0,1) if not given)
% h_local is 3xM (M: number of half vectors)
% density is 1xM
% N (if given) is 3x1. N is assumed to be unit length. Multiple normal is 
% not supported
% There doesn't seem to be a close form of normalization constant for this:
% the result I got from Mathematica is 
%   C = 1 / (2*pi*sqrt(2)*sigma*F(sigma/sqrt(2)))
% where F is "Dawson integral" that involves Erf function. When sigma is
% 0.1, this ~ 1/0.062623
function density = pGaussian(h_local, sigma, N)
    % Compute normalization constant
    C = 1 / (2*pi*sqrt(2)*sigma*mfun('dawson', sigma/sqrt(2)));
    
    % Implemented using
    % http://renderman.pixar.com/view/cook-torrance-shader
    if (isempty(N))
        n_dot_h = h_local(3,:);
    else
        n_dot_h = sum(h_local .* repmat(N, [1 size(h_local, 2)]));
    end
    alpha = acos(n_dot_h);
    density = C * exp(-alpha.^2/(2*sigma^2));
end