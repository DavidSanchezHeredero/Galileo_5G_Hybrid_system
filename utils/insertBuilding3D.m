function insertBuilding3D(map3D, xb, yb, zb, zt, xmax, xmin, ymax, ymin)
% insertBuilding3D - Insert a 3D building footprint with varying roof height into occupancyMap3D
%   map3D : occupancyMap3D object
%   xb, yb : XY coordinates (base and top)
%   zb, zt : Z values of base and top (can vary)
%   xmax, xmin, ymax, ymin : ENU limits

    % Ensure inputs are columns
    xb = xb(:); yb = yb(:); zb = zb(:); zt = zt(:);

    % Map resolution
    res = 1 / map3D.Resolution; % meters per voxel

    % Grid XY
    xVec = xmin:res:xmax;
    yVec = ymin:res:ymax;
    [Xg, Yg] = meshgrid(xVec, yVec);
    pts2D = [Xg(:), Yg(:)];

    % Footprint mask
    inFootprint = inpolygon(pts2D(:,1), pts2D(:,2), xb, yb);
    xyInside = pts2D(inFootprint,:);

    if isempty(xyInside)
        return; % No overlap
    end

    % Use MINIMUM zb as common ground level
    groundZ = min(zb);

    % Interpolate roof height (zt) at every (X,Y)
    % Remove duplicate (xb,yb) vertices by averaging roof height when needed.
    [xyUnique, ~, ic] = unique([xb yb], 'rows', 'stable');
    ztMean = accumarray(ic, zt, [], @mean);
    
    xb = xyUnique(:,1);
    yb = xyUnique(:,2);
    zt = ztMean;

    Froof = scatteredInterpolant(xb, yb, zt, 'linear', 'nearest');
    roofHeights = Froof(xyInside(:,1), xyInside(:,2)); % roof height at each (X,Y)

    % For each (X,Y), generate Z levels up to roof
    allXYZ = [];

    for i = 1:size(xyInside,1)
        xi = xyInside(i,1);
        yi = xyInside(i,2);
        ziTop = roofHeights(i);

        % Vertical range from groundZ to ziTop
        zLevels = groundZ:res:ziTop;
        if isempty(zLevels) % handle thin roofs
            zLevels = groundZ;
        end

        % Create all voxels for this (x,y)
        Xi = xi * ones(numel(zLevels),1);
        Yi = yi * ones(numel(zLevels),1);
        Zi = zLevels(:);

        allXYZ = [allXYZ; [Xi Yi Zi]]; %#ok<AGROW>
    end

    % Insert only into free voxels
    occ = getOccupancy(map3D, allXYZ);
    freeIdx = occ < map3D.OccupiedThreshold;

    if any(freeIdx)
        setOccupancy(map3D, allXYZ(freeIdx,:), 1);
    end
end
