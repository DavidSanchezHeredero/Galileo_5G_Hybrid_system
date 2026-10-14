%% projPos.m
% =========================================================================
% Author:  David Sanchez-Heredero Martinez
% Date:    2025-09-09
% Project: GNSS-5G Urban Hybrid Simulation
%
% 5G NR Base Station Mapper Library
%
% -------------------------------------------------------------------------
% DESCRIPTION:
%   Get satellite ENU projection into map bounding "slab" intersection
%
% INPUT:
%   xmin : minimum value for x axis - [m]
%   xmax : maximum value for x axis - [m]
%   ymin : minimum value for y axis - [m]
%   ymax : maximum value for y axis - [m]
%   zmin : minimum value for z axis - [m]
%   zmax : maximum value for z axis - [m]
%   dir    : [x,y,z] dimension vector. Represents the line draw by the satellite and the receiver
%   rxENU  : [x,y,z] dimension vector. Represents the receiver position
%
% OUTPUT:
%   xt: projected x position.
%   yt: projected x position
%   zt: projected z position
%   fit: [True/False] - boolean indicating the function success
%
% USAGE EXAMPLE:
%   % Get the projected position of a sv into the 3Dmap.
%   [xt,yt,zt,fit] = projPos(xmin, xmax, ymin, ymax, zmin, zmax, dir,rxENU)
%
% SEE ALSO:
%   galileoSimulation.m - environmentSetup.m
% =========================================================================
function [xt,yt,zt,fit] = projPos(xmin, xmax, ymin, ymax, zmin, zmax, dir,rxENU)
    x = dir(1); y = dir(2); z = dir(3);
    % Find the minimum t parameter.
    t = zeros(5,1);
    t(1) = (zmax - rxENU(3)) / z;
    t(2) = (xmin - rxENU(1)) / x;
    t(3) = (xmax - rxENU(1)) / x;
    t(4) = (ymin - rxENU(2)) / y;
    t(5) = (ymax - rxENU(2)) / y;
    [value,index]   = min(abs(t));
    value           = t(index);
    xt=(x)*value+rxENU(1); yt=(y)*value+rxENU(2); zt=(z)*value+rxENU(3);

    if (xt <= xmax && xt >= xmin) && (yt <= ymax && yt >= ymin) && (zt <= zmax && zt >= zmin)
        fit=true;
    else
        fit=false;
    end
end
