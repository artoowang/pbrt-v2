% h_local is in tangent space (i.e., shading normal N is (0,0,1))
% h_local is 3xM (M: number of half vectors)
% density is 1xM
% There doesn't seem to be a close form of normalization constant for this:
% the result I got from Mathematica is 
%   C = 1 / (2*pi*sqrt(2)*sigma*F(sigma/sqrt(2)))
% where F is "Dawson integral" that involves Erf function. When sigma is
% 0.1, this ~ 1/0.062623
function density = pGaussian(h_local, sigma)
    % Compute normalization constant
    C = 1 / (2*pi*sqrt(2)*sigma*mfun('dawson', sigma/sqrt(2)));
    
    % Implemented using
    % http://renderman.pixar.com/view/cook-torrance-shader
    n_dot_h = h_local(3,:);
    alpha = acos(n_dot_h);
    density = C * exp(-alpha.^2/(2*sigma^2));
end