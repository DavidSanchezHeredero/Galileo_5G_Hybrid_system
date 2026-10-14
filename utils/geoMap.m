function geoMap(LoSbStationsENU, groups)
    % geoMap - Plot grouped base stations on a geographic map.
    %
    % Inputs:
    %   LoSbStationsENU - Nx3 ENU base station coordinates.
    %   groups - Cell array whose elements contain the base station indices
    %            assigned to each group.

    % Load reference configuration.
    addpath('./config');  % Make sure MATLAB can find config.m
    run('config');        % Executes config.m and loads the variables into workspace

    fontAxis   = 15;
    fontTitle  = 17;
    fontLegend = 13;
    gridColor  = [0.82 0.82 0.82];

    % Create the geographic figure and axes.
    fig = createThesisFigure('map');
    gx = geoaxes;
    geobasemap(gx, 'satellite');
    hold(gx, 'on');
    title(gx, 'Base Stations Grouped');
    geolimits(gx, [latmin, latmax + 0.005], [lonmin, lonmax]);

    % Define one color per group.
    numGroups = numel(groups);
    colors = lines(numGroups);  % Distinct group colors.

    % Initialize legend entries.
    legendEntries = cell(1, numGroups);

    % Plot each group with its assigned color.
    for k = 1:numGroups
        groupIndices = groups{k};
        groupColor = colors(k, :);
        latitudes = zeros(length(groupIndices), 1);
        longitudes = zeros(length(groupIndices), 1);
        for j = 1:length(groupIndices)
            idx = groupIndices(j);
            enu = LoSbStationsENU(idx, :);
            [lat, lon, ~] = enu2geodetic(enu(1), enu(2), enu(3), refLat, refLon, h0, wgs84Ellipsoid('meters'));
            latitudes(j) = lat;
            longitudes(j) = lon;
        end
        % Plot the group.
        geoscatter(gx, latitudes, longitudes, 36, groupColor, 'filled', 'DisplayName', sprintf('Group %d', k));
        legendEntries{k} = sprintf('Group %d', k);
    end

    % Add legend.
    legend(gx, 'show', 'Location', 'bestoutside');
    styleThesisFigure(fig, fontAxis, fontTitle, fontLegend, gridColor, 'map');
end
