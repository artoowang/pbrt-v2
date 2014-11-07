% The circle is radially symmetric: vals defined its value from the center
% to the boundary, uniformly spaced. E.g., vals = [1 2 3 4] creates 4
% circles with values 1, 2, 3, 4, and radii 0.25, 0.5, 0.75, 1
function drawCircle(vals)
    vals = vals(:); % Nx1
    N = length(vals);
    colors = repmat(vals, [1 3]); % Nx3
    radii = linspace(0, 1, N+1);
    radii = radii(2:end); % 1xN
    radiusMultiplier = 100;
    areas = pi * (radii * radiusMultiplier).^2; % 1xN
    scatter(zeros(1, N), zeros(1, N), fliplr(areas), flipud(colors), 'filled');
end