% =========================================================================
% Script: generateReceiverPath.m 
% Author: David Sanchez-Heredero Martinez
% Date: 2025-04-28
% Project: GNSS-5G Urban Hybrid Simulation
% Description:
%   Generates receiver positions for a walking user inside the simulated
%   environment.
%
% Usage:
%   Run environmentSetup.m first with the desired configuration.
% =========================================================================

clc; clear; close all;

%% Configuration
% Initialize ellipsoid and geoid model
addpath('./utils');   % Add utility functions to the path
addpath('./config');  % Make sure MATLAB can find config.m
run('config');        % Executes config.m and loads the variables into workspace


%% 0) Define path waypoints
locationsLLA = [
    40.420022, -3.696504;
    40.420056, -3.695643;
    40.421869, -3.695671
];
receiverHeight = 1.5;    % Receiver height above terrain [m]

% Parse control points to ENU
numWaypoints = size(locationsLLA, 1);
waypointsENU = zeros(numWaypoints, 3);
for i = 1:numWaypoints
    lat = locationsLLA(i,1);
    lon = locationsLLA(i,2);
    H_ground = interpGeolocated(Z, R, lat, lon);
    N = egm96geoid(lat, lon);
    h_ground = H_ground + N;
    [x, y, z] = geodetic2enu(lat, lon, h_ground, refLat, refLon, h0, ellipsoid);
    waypointsENU(i,:) = [x, y, z + receiverHeight];
end

%% 1) Define movement parameters
speed = 1.5;             % Receiver speed [m/s]
samplingPeriod = 0.5;    % Sampling period [s]
samplingDistance = speed * samplingPeriod;

%% 2) Generate trajectory
receiverPositions = [];
for i = 1:(numWaypoints - 1)
    startPoint = waypointsENU(i,:);
    endPoint = waypointsENU(i+1,:);
    segmentVector = endPoint - startPoint;
    segmentDistance = norm(segmentVector);
    numSamples = ceil(segmentDistance / samplingDistance);
    t = linspace(0, 1, numSamples)';
    segmentPoints = startPoint + t .* segmentVector;
    receiverPositions = [receiverPositions; segmentPoints];
end

%% 3) Add jitter
% Noise - [meters]
jitterStd = 0.2;  % 20 cm standard deviation
receiverPositions = receiverPositions + jitterStd * randn(size(receiverPositions));

%% 4) Estimate velocity
% Position differences
deltaPos = diff(receiverPositions);
% Divide by sampling period
receiverVelocity = deltaPos / samplingPeriod;
% Preserve vector length by prepending a zero velocity sample.
receiverVelocity = [zeros(1,3); receiverVelocity];

%% 5) Plot
[latPoly, lonPoly, altPoly]  = enu2geodetic(...
    receiverPositions(:,1),...
    receiverPositions(:,2),...
    receiverPositions(:,3),...
    refLat, refLon, h0, ellipsoid);

figure;
gx = geoaxes;
geobasemap(gx,'satellite');
hold(gx, 'on'); % Hold the geoaxes rather than the figure.
alpha(gx, 0.7);
title(gx, 'Receiver route');
latmin = bbox(1); latmax = bbox(2);
lonmin = bbox(3); lonmax = bbox(4); 
geolimits(gx, [latmin latmax+0.005], [lonmin lonmax]);
geoscatter(gx, latPoly(:), lonPoly(:), 36, 'r', 'filled');

%% 6) Save
samplingRate = 1 / samplingPeriod; % Sampling rate [Hz]
save(rxENUFile, 'receiverPositions', 'receiverVelocity', 'samplingRate');
save(rxLLAFile, 'latPoly', 'lonPoly', 'altPoly', 'samplingRate');
