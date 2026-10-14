%INTERPGEOLOCATED Interpolate value(s) from a georeferenced raster at specific lat-lon.
%
%   h = interpGeolocated(Z, R, lat, lon) returns the value(s) from the 
%   georeferenced grid Z at the specified latitude and longitude coordinates.
%
%   INPUTS:
%       Z   - 2D matrix of geospatial values (e.g., elevation).
%       R   - GeographicCellsReference object referencing Z.
%       lat - Latitude(s) of query point(s) (in degrees).
%       lon - Longitude(s) of query point(s) (in degrees).
%
%   OUTPUT:
%       h   - Value(s) at the specified geographic coordinates. If lat and lon 
%             are arrays, h will match their size.
%
%   EXAMPLE:
%       [Z, R] = readgeoraster('elevation_madrid.tif');
%       lat = 40.4168;
%       lon = -3.7038;
%       h = interpGeolocated(Z, R, lat, lon);
%
%   NOTE:
%       This function uses nearest-neighbor lookup (discrete sampling).
%       Ensure that the input coordinates lie within the bounds of R.
%
%   See also readgeoraster, geographicToDiscrete, GeographicCellsReference
function h = interpGeolocated(Z, R, lat, lon)
  % R is a GeographicCellsReference
  [row, col] = geographicToDiscrete(R, lat, lon);
  h = Z(row, col);
end