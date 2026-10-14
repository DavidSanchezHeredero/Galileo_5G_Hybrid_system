%% config.m

%% File Paths
% Input
bsFile    = fullfile('./data', 'madrid_bs.json');        % JSON with opencellid information
osmFile   = fullfile('./data', 'madrid.osm');            % OSM file with building footprints
demFile   = fullfile('./data', 'coverage.tif');          % DEM GeoTIFF file
galFile   = fullfile('./data', 'galileo.tle');           % Galileo TLE file

% Generated
mapFile     = fullfile('./data/environment', 'map3D_madrid.mat');          % Map of simulation
rxENUFile   = fullfile('./data/environment', 'receiverTrajectory.mat');    % Receiver ENU data
rxLLAFile   = fullfile('./data/environment', 'receiverTrajectoryLLA.mat'); % Receiver LLA data

bsENUFile   = fullfile('./data/environment', 'baseStationsENU.mat');       % Base station ENU file
nrLoSFile   = fullfile('./data/environment', 'dynamicLoSResults.mat');     % NR LoS information
nrSimFile   = fullfile('./data/estimations', 'nr_estimation.mat');         % NR positioning results
nrObsFile   = fullfile('./data/observables', 'nrTDOA.mat');


galLoSFile  = fullfile('./data/observables', 'galileoLoSResults.mat');    % Galileo LoS evaluation
galSimFile  = fullfile('./data/estimations', 'galileo_estimation.mat');   % Galileo PNT results

ekfSimFile  = fullfile('./data/estimations', 'ekalmanfilter_estimation.mat');   % EKF Positioning results


%% ENU initialization
refLat = 40.4200;                             % Reference latitude (Gran Via)
refLon = -3.6980;                             % Reference longitude
latmax = 40.4224;
latmin = refLat - 0.005;
lonmax = -3.6920;
lonmin = -3.69696;
bbox   = [latmin, latmax, lonmin, lonmax];    % Limits of the simulation
% Initialize ellipsoid and geoid model
global ellipsoid;
ellipsoid = referenceEllipsoid('wgs84','m');
[Z, R] = readgeoraster(demFile, OutputType='double');
% Extract orthometric height at reference point - (local function)
H0 = interpGeolocated(Z, R, refLat, refLon);
% Geoid undulation at reference point (EGM96)
N0 = egm96geoid(refLat, refLon);
% Ellipsoidal height of reference (ENU origin)
h0 = H0 + N0;
refLocation = [refLat, refLon, h0];

%% Time definition - begin of the simulation
startTime = datetime(2025,10,16,10,0,0); % 16th of october, 2025, 10:00:00
