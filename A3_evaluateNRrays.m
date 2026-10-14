% =========================================================================
% Script: evaluateDynamicLoS.m
% Description:
%   Evaluates dynamic LoS visibility and ranging while the receiver moves
%   along the simulated trajectory.
% =========================================================================

clc; clear; close all;

%% 0) Configuration
% Initialize ellipsoid and geoid model
addpath('./utils');   % Add utility functions to the path
addpath('./config');  % Make sure MATLAB can find config.m
run('config');        % Executes config.m and loads the variables into workspace

%% 1) Load map, NR stations and receiver path
load(mapFile, 'map3D', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax');
load(bsENUFile, 'baseStationsENU');
load(rxENUFile, 'receiverPositions', 'samplingRate');

%% 2) Initialize results
nBS = size(baseStationsENU,1);
nRX = size(receiverPositions,1);

LoSMatrix = false(nRX, nBS);       % [nRX x nBS] -> true when LoS is available
RangingMatrix = nan(nRX, nBS);     % [nRX x nBS] -> range when LoS is available

%% 3) Evaluate LoS and ranging
for iRX = 1:nRX
    rx = receiverPositions(iRX,:); % Current receiver position

    for iBS = 1:nBS
        tx = baseStationsENU(iBS,:); % Base station position

        % Compute ray direction and range.
        direction = (rx - tx) / norm(rx - tx);
        range = norm(rx - tx);

        % Define sensor pose as position plus identity quaternion.
        sensorPose = [tx, 1, 0, 0, 0];

        % Trace the ray through the occupancy map.
        [intersection, isOccupied] = rayIntersection(map3D, sensorPose, direction, range, true);

        % Determine whether the receiver is visible from the base station.
        if isOccupied && norm(intersection - rx) < 1e-3
            % The ray reached the receiver without obstacles.
            LoSMatrix(iRX,iBS) = true;
            RangingMatrix(iRX,iBS) = range;
        elseif ~isOccupied
            % The ray found no obstacles before the maximum range.
            LoSMatrix(iRX,iBS) = true;
            RangingMatrix(iRX,iBS) = range;
        else
            % The ray was blocked by an obstacle.
            LoSMatrix(iRX,iBS) = false;
            RangingMatrix(iRX,iBS) = NaN;
        end
    end
end

% Remove base stations that are never visible during the trajectory.
baseStationsNeverLoS    = all(~LoSMatrix, 1);  % Logical vector of size [1 x nBS]
indicesLoS              = find(~baseStationsNeverLoS);  % Base stations with at least one LoS epoch

% Keep only base stations with at least one LoS epoch.
LoSMatrix       = LoSMatrix(:, indicesLoS);
RangingMatrix   = RangingMatrix(:, indicesLoS);
LoSbStationsENU = baseStationsENU(indicesLoS,:);

%% 4) Save results
save(nrLoSFile, 'LoSMatrix', 'RangingMatrix', 'receiverPositions', 'samplingRate', 'LoSbStationsENU');

%% 5) Quick visualization
samplingPeriod  = 1/samplingRate;
timeVector      = (0:nRX-1) * samplingPeriod; % Time vector [s]
figure;
imagesc(timeVector, 1:nBS, LoSMatrix');  % Transpose to place base stations on the y-axis
colormap(gray);
xlabel('Time (s)');
ylabel('Base Station Index');
title('LoS Matrix: Epoch vs Base Station');
colorbar;
set(gca, 'YDir', 'normal');  % Keep base station index 1 at the bottom.

% Count visible base stations per epoch.
numLoS = sum(LoSMatrix, 2); % Column vector [nRX x 1]

% Plot the number of visible base stations over time.
figure;
bar(timeVector, numLoS, 'BarWidth', 1, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('Number of base stations with LoS');
title('LoS conditions over time');
grid on;
