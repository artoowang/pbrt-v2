function [integral] = TEST_BECKMANN(alpha, theta_o)
% view vector
V = [sin(theta_o) 0 cos(theta_o)];
% masking (rational approximation for Lambda)
a = 1 / (alpha * tan(theta_o));
if a < 1.6
Lambda = (1 - 1.259*a + 0.396*a^2) / (3.535*a + 2.181*a^2);
else
Lambda = 0;
end
%G = 1 / (1 + Lambda);
G = G1_Beckmann(V', [0 0 1]', alpha);
integral = 0;
dtheta = 0.05;
dphi = 0.05;
for theta = 0:dtheta:pi
for phi = 0:dphi:2*pi
% reflected vector
L = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
% half vector
H = (V + L) / norm(V + L);
% Beckmann distribution
if H(3) > 0
    % angle associated with H
    theta_h = acos(H(3));
    %D = exp(-(tan(theta_h)/alpha)^2) / (pi * alpha^2 * H(3)^4);
    D = D_Beckmann(H', alpha);
else
    continue;
end
%G = G1_Beckmann(V', H', alpha);
% integrate
integral = integral + sin(theta) * D * G / abs(4 * V(3));
end
end
% display integral (should be 1)
integral = integral * dphi * dtheta;
end