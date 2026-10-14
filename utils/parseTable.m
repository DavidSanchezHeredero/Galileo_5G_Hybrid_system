%% parseTable.m
% =========================================================================
% Author: David Sanchez-Heredero Martinez
% Date: 2025-04-23
% Project: GNSS-5G Urban Hybrid Simulation
%
% Building Footprint Parser
% Parses WKT POLYGON / MULTIPOLYGON strings and extracts the outer ring
% longitude-latitude vertex coordinates of a building footprint.
%
% Functions:
% [latPoly, lonPoly] = parseTable(shp)
% - Parse a WKT POLYGON or MULTIPOLYGON string.
% - Return latitude and longitude vectors for the outer polygon ring.
%
% Input Arguments:
% shp     - Character vector or string scalar containing a WKT POLYGON or
%           MULTIPOLYGON description of a building footprint.
%
% Output Arguments:
% latPoly - Column vector of vertex latitudes of the outer ring.
% lonPoly - Column vector of vertex longitudes of the outer ring.
%
% Notes:
% - Only the exterior ring of the first polygon is considered.
% - Coordinates in the WKT are assumed to be given as "lon lat" pairs.
%
% Usage Example:
% % Extract footprint vertices from a WKT polygon string
% wkt = "POLYGON((-3.70379 40.41678,-3.70370 40.41680,-3.70360 40.41675))";
% [latPoly, lonPoly] = parseTable(wkt);
% =========================================================================
function [latPoly, lonPoly] = parseTable(shp)
    % Erase prefix/suffix text and extract only the exterior ring.
    % Remove the "POLYGON((" prefix and matching suffix.
    str = split(shp, "(");
    if contains(str(1),"MULTIPOLYGON")
        k=4;
    else
        k   = 3; % corresponds to the index related to the building premises
    end
    
    % Each face of the building is presented here
    vertices = split(str(k), ",");
    if vertices(end) == ""
        vertices = vertices(1:end-1);
    end
    vertices(end) = erase( vertices(end), ")" );
    
    coords = zeros(numel(vertices), 2);  % Lat Lon each row
    for v = 1:numel(vertices)
        coords(v,:) = str2double(split(vertices(v(:))," "));
    end
    latPoly = coords(:,2);
    lonPoly = coords(:,1);
end
