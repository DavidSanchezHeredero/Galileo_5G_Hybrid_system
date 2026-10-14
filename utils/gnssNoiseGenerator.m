%% gnssNoiseGenerator.m
% =========================================================================
% Author:  David Sanchez-Heredero Martinez
% Date:    2025-09-09
% Project: GNSS-5G Urban Hybrid Simulation
%
% 5G NR Base Station Mapper Library
% Utility function to inspect MATLAB product dependencies.
%
% [noisyRange, pathLoss_dB] = gnssNoiseGenerator(trueRange, visibility)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   Adds Gaussian noise to true range
%
% INPUT:
%   trueRange : real pseudorange from a satellite and a receiver. Expressed
%               in [m]
%   visibility : Indicates the LoS or NLoS conditions.
%                   the standard deviation is 1m for LoS and 5m otherwise
%                   (1=LoS, 0=NLoS, -1=blocked)
%
% OUTPUT:
%   noisyRange: the new pseudorange with the noise effects.
%   pathLoss_dB: path loss in dB for the pseudorange. [NOT IMPLEMENTED]
%
% USAGE EXAMPLE:
%   % Get a simulated pseudorange from the real value.
%   [noisyRange, ~] = gnssNoiseGenerator(trueRange, visibility)
%
% SEE ALSO:
%   galileoSimulation.m
% =========================================================================
function [noisyRange, pathLoss_dB] = gnssNoiseGenerator(trueRange, visibility)
    % GNSSNOISEGENERATOR: Add Gaussian noise to true range.
    % (Path-loss is not used for GNSS code-phase, so return 0)
    sigma_LoS  = 4.6;    % [m]
    sigma_NLoS = 6.5;    % [m]
    switch visibility
        case 1
            noise = sigma_LoS * randn();
        case 0
            noise = sigma_NLoS * randn();
        otherwise
            noise = NaN;
    end
    noisyRange = trueRange + noise;
    pathLoss_dB = 0;
end
