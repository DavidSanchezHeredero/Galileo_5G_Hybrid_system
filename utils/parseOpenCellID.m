%% parseOpenCellID.m
% =========================================================================
% Author:              David Sanchez-Heredero Martinez
% Date:                2025-04-23
% Project:             GNSS-5G Urban Hybrid Simulation
%
% 5G NR Base Station Mapper Library
%   Provides functions to parse OpenCellID GeoJSON, store base station data,
%   and compute ENU positions and heights relative to urban geometry.
%
% Functions:
%   parseOpenCellID(jsonFile, outMat)
%       - Read a GeoJSON from OpenCellID and save to a .mat structure.
%
% Usage Example:
%   % Parse GeoJSON into MAT
%   parseOpenCellID('openCellID.geojson','baseStations.mat');
% =========================================================================

function parseOpenCellID(jsonFile, outMat)
    % parseOpenCellID  Decode OpenCellID GeoJSON and save baseStations.mat
    %   jsonFile : path to GeoJSON
    %   outMat   : output MAT filename (e.g. 'baseStations.mat')

    txt = fileread(jsonFile);
    data = jsondecode(txt);

    N = numel(data.features);
    Lat = zeros(N,1);
    Lon = zeros(N,1);
    Radio = cell(N,1);
    MCC = zeros(N,1);
    NET = zeros(N,1);
    Area = zeros(N,1);
    Samples = zeros(N,1);
    Range = zeros(N,1);
    Updated = zeros(N,1);

    for k = 1:N
        feat = data.features(k);
        coords = feat.geometry.coordinates;
        Lon(k) = coords(1);
        Lat(k) = coords(2);

        props = feat.properties;
        Radio{k}   = props.radio;
        MCC(k)     = props.mcc;
        NET(k)     = props.net;
        Area(k)    = props.area;
        Samples(k) = props.samples;
        Range(k)   = props.range;       % coverage radius, km or m
        Updated(k) = props.updated;     % posix timestamp
    end

    UpdatedTime = datetime(Updated,'ConvertFrom','posixtime','TimeZone','UTC');

    baseStations = struct();
    baseStations.Lat     = Lat;
    baseStations.Lon     = Lon;
    baseStations.Radio   = Radio;
    baseStations.MCC     = MCC;
    baseStations.NET     = NET;
    baseStations.Area    = Area;
    baseStations.Samples = Samples;
    baseStations.Range   = Range;
    baseStations.Updated = UpdatedTime;

    save(outMat,'baseStations');
    fprintf('Saved %d base stations to %s\n', N, outMat);
end  
