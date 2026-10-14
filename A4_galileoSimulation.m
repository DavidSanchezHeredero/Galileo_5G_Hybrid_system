% =========================================================================
% Script: galileoSimulation.m
% Author: David Sanchez-Heredero Martinez
% Date: 2025-04-28
% Project: GNSS-5G Urban Hybrid Simulation
% Description:
%       GNSS - Galileo Simulation with Newton-Based Position Estimation.
% 0) Measurement Parameters
% 1) Build Satellite Scenario
% 2) Pre-allocate Storage
% 3) Loop Over GNSS Epochs
%   3.1) Receiver True Position in ENU
%   3.2) Create Temporary Ground Station at Receiver LLA
%   3.3) Determine Visible Satellites via access()
%   3.4) Extract ECEF-to-ENU for visible satellites
%   3.5) For Each Visible Satellite, Compute Pseudorange 
%       3.5.1) True geometric range
%       3.5.2) Check ray-slab intersection with map3D
%   3.6) Estimate Receiver Position (ECEF)
% 4) Save results
%
% Usage:
%   Run environmentSetup.m and generateReceiverPath.m first.
% =========================================================================

clc; clear; close all;

%% Configuration and Data Load
addpath('./utils');   % Add utility functions to the path
addpath('./config');  % Make sure MATLAB can find config.m
run('config');        % Executes config.m and loads the variables into workspace

demFile   = 'data/coverage.tif';          % DEM GeoTIFF file
load(mapFile, 'map3D', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax');
load(rxENUFile, 'receiverPositions', 'samplingRate');
load(rxLLAFile, 'latPoly', 'lonPoly', 'altPoly');

%% 0) Measurement Parameters      
% User
epochs      = size(receiverPositions,1);
Fs_user     = samplingRate;
% Galileo
Fs_gnss     = Fs_user/2;                    % 1 Hz GNSS update rate
nGnss       = floor(Fs_gnss * (epochs / Fs_user));
% Fixed start time
tGnss = startTime + seconds((0:nGnss-1)'/Fs_gnss);

%% 1) Build Satellite Scenario
startTime = tGnss(1);
stopTime  = tGnss(end);
sc        = satelliteScenario(startTime, stopTime, Fs_gnss);
sats      = satellite(sc, galFile, "OrbitPropagator", "sgp4");
play(sc);
% Retrieve ECEF state history: 3 x nGnss x numSats.
satGeoAll = states(sats, "CoordinateFrame", "ecef");
numSats   = numel(sats);

%% 2) Pre-allocate Storage
svPositions_ENU  = zeros(numSats, 3, nGnss);    % ENU positions of each SV each epoch
svPositions_ECEF = zeros(numSats, 3, nGnss);    % ECEF positions for Newton
trueRange        = nan(numSats, nGnss);         % geometric range
visibilityMap    = zeros(numSats, nGnss);       % 1=LoS, 0=elevation-nlos, -1=blocked
measPseudorange  = nan(numSats, nGnss);         % noisy pseudorange measurements
estENU           = nan(nGnss, 3);                % estimated ENU positions per epoch
estError         = nan(nGnss, 1);                % estimation error per epoch
TsGalileo_pos    = 1 / Fs_gnss;                  % GNSS update period [s]

%% 3) Loop Over GNSS Epochs
for k = 1:nGnss
    fprintf('--- Epoch %d / %d ---\n', k, nGnss);
    t = tGnss(k);

    % 3.1) Receiver True Position in ENU
    idxRx = round(k * Fs_user / Fs_gnss);
    idxRx = max(min(idxRx, epochs), 1);
    rxENU = receiverPositions(idxRx, :);            % ENU [x, y, z]
    rxLLA = [latPoly(idxRx), lonPoly(idxRx), altPoly(idxRx)];

    % 3.2) Create Temporary Ground Station at Receiver LLA
    gs = groundStation(sc, ...
        'Latitude',         rxLLA(1), ...
        'Longitude',        rxLLA(2), ...
        'Altitude',         rxLLA(3), ...
        'MinElevationAngle', 10);

    % 3.3) Determine Visible Satellites via access()
    ac     = access(sats, gs);       % Logical vector numSats x 1
    intvls = accessIntervals(ac);    % table of intervals
    satNames = {sats.Name}';         % 1 x 1 cell
    satNames = satNames{1};          % 1 x numSats string array
    srcNames = intvls.Source;         % K x 1 string array
    [tf, idxFromTable] = ismember(srcNames, satNames);
    visibleSatIdx = unique(idxFromTable(tf));  % Indices of satellites in view

    if isempty(visibleSatIdx)
        warning('Epoch %d: No satellites above elevation mask. Skipping.', k);
        continue;
    end

    %% 3.4) Extract ECEF-to-ENU for visible satellites
    % satGeoAll is 3 x nGnss x numSats, so all satellites at epoch k are
    % recovered as squeeze(satGeoAll(:,k,:)) -> 3 x numSats.
    satECEF_k = squeeze(satGeoAll(:, k, :))';   % numSats x 3
    % Convert all satellites to ENU (will index later)
    [xAll, yAll, zAll] = ecef2enu( ...
        satECEF_k(:,1), satECEF_k(:,2), satECEF_k(:,3), ...
        refLat, refLon, h0, ellipsoid);
    satENU_all = [xAll, yAll, zAll];             % numSats x 3

    %% 3.5) Compute pseudorange for each visible satellite
    visibleCount    = numel(visibleSatIdx);
    satPos_ECEF_vis = zeros(visibleCount, 3);
    rho_meas        = nan(visibleCount, 1);
    trans = zeros(visibleCount,3);
    for j = 1:visibleCount
        s       = visibleSatIdx(j);
        rtxENU  = satENU_all(s, :);        % ENU sat position
        svPositions_ENU(s, :, k)    = rtxENU; 

        svPositions_ECEF(s, :, k)   = satECEF_k(s, :);
        % 3.5.1) True geometric range
        rTrue = norm(rtxENU - rxENU);
        trueRange(s, k) = rTrue;

        % 3.5.2) Check ray-slab intersection with map3D
        % Project satellite ENU to map boundary if needed
        direction = (rxENU - rtxENU) / norm(rxENU - rtxENU);
        [xt, yt, zt, fit] = projPos(xmin, xmax, ymin, ymax, zmin, zmax, direction,rxENU);
        Pos = [xt, yt, zt];
        if ~fit
            fprintf(   'Could not fit sat %d within the map3d box.',s);
        end
      
        txENU = [xt, yt, zt];
        trans(s,:) = [xt, yt, zt];
        maxRange  = norm(rxENU - txENU);
        sensorPose = [txENU, 1, 0, 0, 0];  % [x,y,z, qw,qx,qy,qz] with identity quaternion
        
        [intersection, isOcc] = rayIntersection(map3D, sensorPose, direction, maxRange, true);
        if isOcc  == -1  || isOcc  == -0
            % Satellite visible (LoS)
            visibilityMap(s, k) = 1;
            % Add GNSS-style noise
            [rNoisy, ~] = gnssNoiseGenerator( ...
                rTrue,          ... % true range [m]
                1);                 % visibility flag
            measPseudorange(s, k) = rNoisy;
            rho_meas(j) = rNoisy;
            satPos_ECEF_vis(j, :) = satECEF_k(s, :);
            fprintf('  Sat %d LoS: true=%.1f m, noisy=%.1f m\n', s, rTrue, rNoisy);
        else
            % Blocked by building
            visibilityMap(s, k) = -1;
            measPseudorange(s, k) = NaN;
            fprintf('  Sat %d blocked by building at epoch %d.\n', s, k);
            continue;
        end
    end

    %% 3.6) Estimate Receiver Position (ECEF)
    indices     = find(visibilityMap(:, k) == 1);
    rho         = measPseudorange(indices, k);

    % Must have at least 4 satellites for a unique solution
    if numel(rho) < 4
        warning('  Epoch %d: insufficient visible satellites (%d). Skipping estimation.', k, numel(rho));
        continue;
    end

    estimation = receiverposition(trueRange(indices,k),satECEF_k(indices,:));

    % p now contains [x_est; y_est; z_est; clk_bias_est]
    [estX, estY, estZ] = deal(estimation(1), estimation(2), estimation(3));
    fprintf('  Estimated ECEF: [%.1f, %.1f, %.1f]\n', ...
            estX, estY, estZ);

    % Convert estimated ECEF to ENU for error computation.
    %[est_e, est_n, est_u] = ecef2enu(estX, estY, estZ, ...
    % Convert LLA to ENU.
    [est_e, est_n, est_u] = geodetic2enu(estX, estY, estZ, ...
                                     refLat, refLon, h0, ellipsoid);
    estENU(k, :) = [est_e, est_n, est_u];
    % Compute estimation error (ENU norm)
    err = norm(estENU(k,:) - rxENU);
    estError(k) = err;
    fprintf('  True ENU:     [%.1f, %.1f, %.1f]\n', rxENU);
    fprintf('  Est. ENU:     [%.1f, %.1f, %.1f], error = %.2f m\n', ...
            est_e, est_n, est_u, err);
end

%% 4) Save results
svPositionsECEF = satGeoAll;
save(galLoSFile, 'visibilityMap', 'svPositionsECEF', 'svPositions_ENU', 'trueRange', 'measPseudorange', 'Fs_gnss');
save(galSimFile, 'estENU', 'estError', 'TsGalileo_pos');
