% =========================================================================
% Author: David Sanchez-Heredero Martinez
% Date: 2025-04-28
% Project: GNSS-5G Urban Hybrid Simulation
% Description:
%   Compute ENU positions for buildings and base stations in Madrid.
%   1) Load OSM building footprints and display.
%   2) Load DEM GeoTIFF and compute ellipsoidal heights.
%   3) Parse building polygons, convert to ENU collision boxes.
%   4) Load base station ENU positions.
%   5) Visualize a 3D building model using wall and roof patches.
% Usage:
%   - Place 'madrid.osm' and 'coverage.tif' in current folder.
%   - Adjust reference location and filenames as needed.
% =========================================================================
clc; clear; close all;

%% Configuration
addpath('./utils');   % Add utility functions to the path
addpath('./config');  % Make sure MATLAB can find config.m
run('config');        % Executes config.m and loads the variables into workspace

%% 0) Map OpenCellID base stations
parseOpenCellID(bsFile, "data/baseStations.mat");

%% 1) Load building footprints
fprintf('Loading data...\n');
buildingsLayer = readgeotable(osmFile, Layer='buildings');
fprintf('  Total number of buildings obtained from the OSM file: %d\n', height(buildingsLayer));

%% 2) Load DEM and reference heights
[Z, R] = readgeoraster(demFile, OutputType='double');
% Extract orthometric height at reference point
H0 = interpGeolocated(Z, R, refLat, refLon);
% Geoid undulation at reference point (EGM96)
N0 = egm96geoid(refLat, refLon);
% Ellipsoidal height of reference (ENU origin)
h0 = H0 + N0;
refLocation = [refLat, refLon, h0];


%% 3) Compute and visualize ENU model for a single building
% Extract geodetic centroid coordinates.
lat = buildingsLayer.Centroid.Latitude;
lon = buildingsLayer.Centroid.Longitude;

% bbox = [latMin latMax lonMin lonMax]
inMask = (lat >= bbox(1) & lat <= bbox(2)) & ...
       (lon >= bbox(3) & lon <= bbox(4));

% Indices of buildings inside the configured bounds.
inIdx  = find(inMask); fprintf('  Buildings inside the area: %d\n', numel(inIdx));
outIdx = find(~inMask);

%% 4) Convert buildings to ENU and create occupancyMap3D
rowGT       = buildingsLayer(:,:);
T           = geotable2table(rowGT, ["Shape"]);

% Create a ground plane in the local ENU frame.
% Map limits for our analysis
[xmax,ymax,~]=geodetic2enu(latmax, lonmax, 700, refLat,refLon,700,ellipsoid);
[xmin,ymin,~]=geodetic2enu(latmin, lonmin, 700, refLat,refLon,700,ellipsoid);

% Voxel resolution [cells/m].
resolution  = 1;
map3D       = occupancyMap3D(resolution);

% Initialize vertical-limit accumulators.
allZb = [];
allZt = [];

fprintf('  Processing buildings...\n');

for j = 1:sum(inMask)
    i = inIdx(j);
    % Extract latitude-longitude vertices.
    [latPoly, lonPoly] = parseTable(T.Shape(i));
    %centroids = geoscatter(gx, latPoly, lonPoly, 36, 'g', 'filled');

    % Convert building vertices to ENU.
    [xb, yb, zb, xt, yt, zt] = building2ENU(...
        latPoly,lonPoly,R,Z,...
        buildingsLayer.MaxHeight(i),...
        refLocation,ellipsoid ...
        );

    % Insert this building prism into the occupancy grid.
    insertBuilding3D(map3D, xb, yb, zb, zt,xmax,xmin,ymax,ymin);

    % Update base and top height samples.
    allZb = [allZb; zb(:)];
    allZt = [allZt; zt(:)];

    % Print progress every 50 buildings.
    if mod(j,50) == 0
        fprintf('    %d/%d buildings processed...\n', j, numel(inIdx));
    end
end
% Derive vertical map limits after all buildings have been processed.
zmin = min(allZb);
zmax = max(allZt);

fprintf('  All buildings where processed.\n\n');

%% 5) BASE STATIONS
fprintf('Loading Base Stations...\n');

if ~exist("data/baseStations.mat","file")
    error('Error: First, execute baseStationMapperLib.m');
end
load('data/baseStations.mat');

lat = baseStations.Lat(:);
lon = baseStations.Lon(:);
    
fprintf('  Total number of Base Stations found: %d\n', numel(lat));

% bbox = [latMin latMax lonMin lonMax]
inMask = (lat >= bbox(1) & lat <= bbox(2)) & ...
       (lon >= bbox(3) & lon <= bbox(4));

% Index vector
inIdx  = find(inMask);

fprintf('  Inside the area where found: %d\n', numel(inIdx));

%% 6) BS TO ENU
%% Terrain heights and geoid per vertex
fprintf('  Obtaining Base Station ENU\n');
% Preallocate matrix for ENU coordinates
baseStationsENU = zeros(numel(inIdx), 3);

% Keep base stations contained within the configured bounds.
lat_bs = lat(inIdx);
lon_bs = lon(inIdx);

% Interpolate heights from the DEM
hGround_bs = arrayfun(@(la, lo) interpGeolocated(Z, R, la, lo), lat_bs, lon_bs);

% Convert geodetic coordinates to ENU with respect to refLocation.
[xbs, ybs, ~] = geodetic2enu(lat_bs, lon_bs, hGround_bs, refLat, refLon, h0, ellipsoid);
% Build geographic query points for the building layer.
queryPoint = geopointshape(lat, lon);
queryPoint.GeographicCRS = buildingsLayer.Shape(1).GeographicCRS;

% Assign a plausible antenna height for each base station.
fprintf('    .');
for k = 1:numel(inIdx)
    % ENU coordinates
    x = xbs(k);
    y = ybs(k);
    z = randi([3 20]);       % Add randomness because the real antenna height is unknown.

    shape = buildingsLayer.Shape(k);
    if isinterior(shape, queryPoint)
        isInside = true;
        z = z + buildingsLayer.MaxHeight(k);
    else
        z = z + 20;
    end

    % Store the ENU base station coordinates.
    baseStationsENU(k, :) = [x, y, z];
    
    fprintf('.');
    % Progress log.
    if mod(j,10) == 0
        fprintf('\n');
    end
end
fprintf('\n');

%% Store map3D, limits of the map and base stations
save(mapFile, 'map3D', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax');
save(bsENUFile, 'baseStationsENU');

toc; % <-- Show real time

fprintf('Script succesfully executed.\n');
