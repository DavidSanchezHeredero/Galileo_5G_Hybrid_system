
% =========================================================================
% Script: galileoSimulation.m
% Project: GNSS-5G Urban Hybrid Simulation
%
% Galileo simulation with:
%   1. Satellite propagation using SGP4
%   2. Elevation-based satellite visibility
%   3. Signal propagation time correction
%   4. Earth rotation correction (Sagnac)
%   5. Urban LoS verification
%   6. Pseudorange generation
%   7. GNSS receiver position estimation
%
% Usage:
%   Run environmentSetup.m and generateReceiverPath.m first.
% =========================================================================

clc; clear; close all;

%% Configuration and Data Load

addpath('./utils');
addpath('./config');

run('config');

load(mapFile, ...
    'map3D', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax');

load(rxENUFile, 'receiverPositions', 'samplingRate');

load(rxLLAFile, 'latPoly', 'lonPoly', 'altPoly');

%% 0) Measurement Parameters

epochs  = size(receiverPositions,1);
Fs_user = samplingRate;

% GNSS update rate
Fs_gnss = Fs_user/2;

% Include initial sample (t = 0)
duration = (epochs-1)/Fs_user;
nGnss = floor(duration*Fs_gnss) + 1;

% GNSS epoch timestamps
tGnss = startTime + seconds((0:nGnss-1)'/Fs_gnss);

% Physical constants
c = physconst('LightSpeed');

omegaEarth = 7.2921151467e-5; % rad/s

minElevation = 10;             % degrees

%% 1) Build Satellite Scenario

startTime = tGnss(1);
stopTime  = tGnss(end);

sc = satelliteScenario( ...
    startTime, stopTime, 1/Fs_gnss);

sats = satellite(sc, galFile, ...
    "OrbitPropagator", "sgp4");

numSats = numel(sats);

% Optional satellite visualization
% play(sc);

% Retrieve satellite positions and velocities.
% Dimensions: 3 x nGnss x numSats
[satGeoAll, satVelAll] = states( ...
    sats, "CoordinateFrame", "ecef");

%% 2) Pre-allocate Storage

svPositions_ENU = nan(numSats,3,nGnss);

% Satellite ECEF positions at reception time
svPositions_ECEF = nan(numSats,3,nGnss);

% Satellite ECEF positions at transmission time
svPositions_TX_ECEF = nan(numSats,3,nGnss);

% Satellite coordinates corrected to reception ECEF frame
svPositions_Sagnac_ECEF = nan(numSats,3,nGnss);

trueRange = nan(numSats,nGnss);

% Visibility:
%  1 = LoS
%  0 = Below elevation mask
% -1 = Blocked by buildings
visibilityMap = zeros(numSats,nGnss);

measPseudorange = nan(numSats,nGnss);

% Signal propagation diagnostics
signalTravelTime = nan(numSats,nGnss);
earthRotation = nan(numSats,nGnss);

% Position estimation
estENU = nan(nGnss,3);
estError = nan(nGnss,1);

TsGalileo_pos = 1/Fs_gnss;

% Maximum urban ray length
mapDiagonal = norm( ...
    [xmax-xmin, ymax-ymin, zmax-zmin]);

%% 3) Loop Over GNSS Epochs

for k = 1:nGnss

    fprintf('--- Epoch %d / %d ---\n', k, nGnss);

    %% 3.1) Receiver True Position

    idxRx = 1 + round((k-1)*Fs_user/Fs_gnss);

    idxRx = min(max(idxRx,1),epochs);

    rxENU = receiverPositions(idxRx,:);

    rxLLA = [ ...
        latPoly(idxRx), ...
        lonPoly(idxRx), ...
        altPoly(idxRx)];

    % Receiver ECEF coordinates at reception time
    [rxX,rxY,rxZ] = geodetic2ecef( ...
        ellipsoid, ...
        rxLLA(1),rxLLA(2),rxLLA(3));

    rxECEF = [rxX,rxY,rxZ];

    %% 3.2) Satellite States at Reception Time

    % Get positions and velocities of every satellite
    satECEF_k = reshape( ...
        satGeoAll(:,k,:),3,numSats).';

    satVel_k = reshape( ...
        satVelAll(:,k,:),3,numSats).';

    %% 3.3) Compute Satellite Positions at Transmission
    %
    % Satellite positions from SGP4 correspond to reception time.
    %
    % For each satellite:
    %   1. Estimate signal travel time
    %   2. Estimate satellite position at transmission
    %   3. Rotate into reception ECEF frame (Sagnac)
    %   4. Recalculate geometric distance
    %
    % Three iterations are used.

    satTX = zeros(numSats,3);
    satCorrected = zeros(numSats,3);

    rhoGeometric = zeros(numSats,1);
    tauAll = zeros(numSats,1);
    thetaAll = zeros(numSats,1);

    for s = 1:numSats

        posRX = satECEF_k(s,:);
        velRX = satVel_k(s,:);

        % Initial geometric distance
        tau = norm(posRX-rxECEF)/c;

        for iter = 1:3

            % Approximate position at signal transmission
            % using first-order orbital propagation.
            %
            % ECEF velocity is expressed in m/s.
            posTX = posRX - velRX*tau;

            % Earth rotation angle during signal travel
            theta = omegaEarth*tau;

            % Transform satellite coordinates from
            % transmission ECEF to reception ECEF
            R = [ cos(theta),  sin(theta), 0;
                 -sin(theta),  cos(theta), 0;
                  0,           0,          1 ];

            posCorrected = (R*posTX.').';

            % Distance accounting for Earth rotation
            rho = norm(posCorrected-rxECEF);

            % Update propagation time
            tau = rho/c;

        end

        % Recalculate final values with converged tau
        posTX = posRX - velRX*tau;
        theta = omegaEarth*tau;

        R = [ cos(theta),  sin(theta), 0;
             -sin(theta),  cos(theta), 0;
              0,           0,          1 ];

        posCorrected = (R*posTX.').';

        rho = norm(posCorrected-rxECEF);

        % Store calculated values
        satTX(s,:) = posTX;
        satCorrected(s,:) = posCorrected;

        rhoGeometric(s) = rho;
        tauAll(s) = rho/c;
        thetaAll(s) = theta;

    end

    % Save diagnostics
    svPositions_ECEF(:,:,k) = satECEF_k;

    svPositions_TX_ECEF(:,:,k) = satTX;

    svPositions_Sagnac_ECEF(:,:,k) = satCorrected;

    signalTravelTime(:,k) = tauAll;

    earthRotation(:,k) = thetaAll;

    %% 3.4) Satellite Elevation and Visibility
    %
    % No groundStation or access objects are created.
    % Elevation is calculated directly for the receiver
    % location at this epoch.
    %
    % Satellite vectors are expressed in the receiver's
    % ECEF frame at reception.

    % Vector from receiver to satellites
    satDiff = satCorrected - rxECEF;

    % Unit vector pointing to local Up
    lat = deg2rad(rxLLA(1));
    lon = deg2rad(rxLLA(2));

    upVector = [ ...
        cos(lat)*cos(lon), ...
        cos(lat)*sin(lon), ...
        sin(lat)];

    % Geometric distances
    distances = vecnorm(satDiff,2,2);

    % Elevation above local horizontal plane
    sinElevation = (satDiff*upVector.') ./ distances;

    % Protect against numerical errors
    sinElevation = max(-1,min(1,sinElevation));

    elevation = asind(sinElevation);

    % Apply elevation mask
    visibleSatIdx = find(elevation >= minElevation);

    if isempty(visibleSatIdx)

        warning('Epoch %d: no satellites above elevation mask.',k);

        continue;

    end

    %% 3.5) Convert Satellite Positions to ENU

    % Use satellite positions corrected to reception ECEF
    [xAll,yAll,zAll] = ecef2enu( ...
        satCorrected(:,1), ...
        satCorrected(:,2), ...
        satCorrected(:,3), ...
        refLat,refLon,h0,ellipsoid);

    satENU_all = [xAll,yAll,zAll];

    %% 3.6) Urban LoS and Pseudorange Generation

    visibleCount = numel(visibleSatIdx);

    fprintf('  Above elevation mask: %d satellites\n', ...
        visibleCount);

    for j = 1:visibleCount

        s = visibleSatIdx(j);

        rtxENU = satENU_all(s,:);

        svPositions_ENU(s,:,k) = rtxENU;

        % Geometric range including signal propagation
        % and Earth rotation
        rTrue = rhoGeometric(s);

        trueRange(s,k) = rTrue;

        % Direction from receiver towards satellite
        direction = rtxENU-rxENU;
        direction = direction/norm(direction);

        % Start the ray at the receiver location
        sensorPose = [rxENU,1,0,0,0];

        % Urban obstruction test.
        % Trace sufficiently far to leave the local map.
        maxRange = min(rTrue,2*mapDiagonal);

        [~,isOcc] = rayIntersection( ...
            map3D, ...
            sensorPose, ...
            direction, ...
            maxRange, ...
            true);

        if isOcc == 0 || isOcc == -1

            % Satellite with direct LoS
            visibilityMap(s,k) = 1;

            % Generate noisy pseudorange
            [rNoisy,~] = gnssNoiseGenerator( ...
                rTrue,1);

            measPseudorange(s,k) = rNoisy;

            fprintf(['  Sat %d LoS: ', ...
                'elev=%.1f deg, ', ...
                'range=%.1f m, ', ...
                'noisy=%.1f m\n'], ...
                s,elevation(s),rTrue,rNoisy);

        else

            % Satellite blocked by buildings
            visibilityMap(s,k) = -1;

            measPseudorange(s,k) = NaN;

            fprintf('  Sat %d blocked at epoch %d.\n',s,k);

        end

    end

    %% 3.7) Estimate Receiver Position

    % Select urban LoS satellites
    indices = find(visibilityMap(:,k) == 1);

    % Use measured pseudoranges, NOT ideal ranges
    rho = measPseudorange(indices,k);

    % IMPORTANT:
    % receiverposition already models Sagnac.
    %
    % Therefore use satellite ECEF positions at
    % transmission WITHOUT manually rotating them.
    satPos = satTX(indices,:);

    if numel(rho) < 4

        warning(['Epoch %d: insufficient LoS ', ...
            'satellites (%d).'],k,numel(rho));

        continue;

    end

    % Receiver position in geodetic coordinates:
    % [latitude, longitude, altitude]
    estLLA = receiverposition(rho,satPos);

    % Convert estimated LLA into local ENU
    [est_e,est_n,est_u] = geodetic2enu( ...
        estLLA(1),estLLA(2),estLLA(3), ...
        refLat,refLon,h0,ellipsoid);

    estENU(k,:) = [est_e,est_n,est_u];

    % Position estimation error
    err = norm(estENU(k,:)-rxENU);

    estError(k) = err;

    fprintf('  True ENU: [%.2f %.2f %.2f]\n',rxENU);

    fprintf('  Est. ENU: [%.2f %.2f %.2f]\n',estENU(k,:));

    fprintf('  Position error: %.3f m\n',err);

end

%% 4) Save Results

% Original satellite state history
svPositionsECEF = satGeoAll;

save(galLoSFile, ...
    'visibilityMap', ...
    'svPositionsECEF', ...
    'svPositions_ENU', ...
    'svPositions_ECEF', ...
    'svPositions_TX_ECEF', ...
    'svPositions_Sagnac_ECEF', ...
    'trueRange', ...
    'measPseudorange', ...
    'signalTravelTime', ...
    'earthRotation', ...
    'Fs_gnss');

save(galSimFile, ...
    'estENU', ...
    'estError', ...
    'TsGalileo_pos');

%% 5) Summary

validErrors = estError(isfinite(estError));

fprintf('\n========== GALILEO SIMULATION ==========\n');

fprintf('Total GNSS epochs: %d\n',nGnss);

fprintf('Successful estimates: %d\n',numel(validErrors));

if ~isempty(validErrors)

    fprintf('Mean error: %.3f m\n',mean(validErrors));

    fprintf('Median error: %.3f m\n',median(validErrors));

    fprintf('RMSE: %.3f m\n',sqrt(mean(validErrors.^2)));

    fprintf('Maximum error: %.3f m\n',max(validErrors));

end

fprintf('========================================\n');
