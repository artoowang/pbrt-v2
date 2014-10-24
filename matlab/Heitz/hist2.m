% xLimit, yLimit: 2-vector specifying the min/max value of x/y
function count = hist2(x, y, xLimit, yLimit, xRes, yRes)
    % Normalize the specified range to [0,1]
    x = (x - xLimit(1)) / (xLimit(2) - xLimit(1));
    y = (y - yLimit(1)) / (yLimit(2) - yLimit(1));
    
    % Compute the index on x/y axes (0-based)
    idx = max(min(floor(x * xRes), xRes-1), 0);
    idy = max(min(floor(y * yRes), yRes-1), 0);
    
    % Compute linear index (1-based)
    ids = idx * yRes + idy + 1;
    
    % Count the numbers in them
    uniqIds = unique(ids);
    idCounts = histc(ids, uniqIds);
    
    % Generate histogram
    count = zeros(yRes, xRes);
    count(uniqIds) = idCounts;
end