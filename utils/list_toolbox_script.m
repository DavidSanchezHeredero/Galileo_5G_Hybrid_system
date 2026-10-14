%% list_toolbox_script.m
% =========================================================================
% Author:  David Sanchez-Heredero Martinez
% Date:    2025-09-09
% Project: GNSS-5G Urban Hybrid Simulation
%
% 5G NR Base Station Mapper Library
% Utility function to inspect MATLAB product dependencies.
%
% list_toolbox_script(fichero)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   Analyses a MATLAB file and prints to the command window the MATLAB
%   products (MATLAB, toolboxes, etc.) required to run it.
%
%   Internally, this function uses matlab.codetools.requiredFilesAndProducts
%   to detect dependencies.
%
% INPUT:
%   fichero : char vector or string scalar
%             Full or relative path to the .m file to be analysed.
%
% OUTPUT:
%   This function does not return any variables. It prints a formatted list
%   of required products/toolboxes to the command window.
%
% USAGE EXAMPLE:
%   % Analyse a script and list required toolboxes
%   list_toolbox_script("B:\code\environmentSetup.m")
%
% SEE ALSO:
%   matlab.codetools.requiredFilesAndProducts, ver
% =========================================================================
function list_toolbox_script(fichero)
    if nargin < 1
        error('Use: list_toolbox_script(''fichero.m'')');
    end
    [~, productos] = matlab.codetools.requiredFilesAndProducts(fichero);
    fprintf('Toolbox needed for %s:\n', fichero);
    for k = 1:numel(productos)
        fprintf(' - %s\n', productos(k).Name);
    end
end
