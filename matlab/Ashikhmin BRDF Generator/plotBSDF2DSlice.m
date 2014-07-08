% bsdfEval: a function of (wi, wr), where wi is 3xN matrix, wr is 3x1
%           vector
function [bsdfVals, theta_is] = plotBSDF2DSlice(bsdfEval, nTheta, wr)
    [~, phi_r] = vector2sph(wr);
    phi_i = mod(phi_r + pi, 2*pi);
    theta_is = linspace(0, pi/2, nTheta);
    wis = sph2vector(theta_is, repmat(phi_i, [1 nTheta]));
    bsdfVals = bsdfEval(wis, wr);
    reshape(bsdfVals, size(theta_is));
    polar(pi/2 - theta_is, bsdfVals);
    ticks = findall(gcf, 'type', 'text');
    delete(ticks);
end