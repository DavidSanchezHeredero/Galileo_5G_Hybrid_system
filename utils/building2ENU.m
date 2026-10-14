%% Convert building footprint to ENU
function [xb, yb, zb, xt, yt, zt] = building2ENU(latPoly,lonPoly,R,Z,height,refLocation,ellipsoid)
    %% Terrain and geoid heights for each vertex
    nV     = numel(latPoly);
    H_poly = nan(nV,1);
    latlim = R.LatitudeLimits;     % [minLat maxLat]
    lonlim = R.LongitudeLimits;    % [minLon maxLon]
    
    lat0    = refLocation(1);
    lon0    = refLocation(2);
    alt0    = refLocation(3);

    for k = 1:nV
        la = latPoly(k);
        lo = lonPoly(k);
        % Check whether the vertex lies inside the DEM coverage.
        if la >= latlim(1) && la <= latlim(2) && lo >= lonlim(1) && lo <= lonlim(2)
            % Interpolate the orthometric height.
            H_poly(k) = interpGeolocated(Z, R, la, lo);
        else
            warning('Vertex %d outside DEM limits: lat=%.6f, lon=%.6f', k, la, lo);
            H_poly(k) = NaN;    % Default value for out-of-bounds vertices.
        end
    end
    H_poly = arrayfun(@(la,lo) interpGeolocated(Z,R,la,lo), latPoly, lonPoly);
    N_poly = egm96geoid(latPoly, lonPoly);
    hBase  = H_poly + N_poly;             % Ellipsoidal ground height at each vertex.
    
    %% Building height
    bldHeight = height;
    hTop      = hBase + bldHeight;        % Ellipsoidal roof height.
    
    %% Convert base and top vertices to ENU
    [xb, yb, zb] = geodetic2enu(latPoly, lonPoly, hBase, ...
                                lat0, lon0, alt0, ellipsoid);
    [xt, yt, zt] = geodetic2enu(latPoly, lonPoly, hTop, ...
                                lat0, lon0, alt0, ellipsoid);
end
