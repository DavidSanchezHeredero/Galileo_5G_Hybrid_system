function figSizeCm = thesisFigureSize(sizeSpec)
%THESISFIGURESIZE Return a thesis-oriented figure size in centimeters.
%   FIGSIZECM = THESISFIGURESIZE(SIZESPEC) accepts either a preset name or a
%   numeric [width height] pair in centimeters.

    if nargin == 0 || isempty(sizeSpec)
        sizeSpec = 'single';
    end

    if isnumeric(sizeSpec)
        validateattributes(sizeSpec, {'numeric'}, {'vector', 'numel', 2, 'positive'});
        figSizeCm = double(reshape(sizeSpec, 1, 2));
        return;
    end

    sizeKey = lower(string(sizeSpec));
    switch sizeKey
        case {"single", "default"}
            figSizeCm = [13.5 9.0];
        case {"wide", "timeseries"}
            figSizeCm = [14.5 9.0];
        case {"map", "geo"}
            figSizeCm = [14.5 10.0];
        case {"stacked", "tall"}
            figSizeCm = [14.0 12.0];
        case {"two-panel", "two_panel"}
            figSizeCm = [15.5 9.0];
        case {"square"}
            figSizeCm = [12.5 12.5];
        case {"grid", "grid-wide", "grid_wide"}
            figSizeCm = [14.5 8.8];
        case {"three-d", "three_d", "3d"}
            figSizeCm = [14.5 10.5];
        otherwise
            error('thesisFigureSize:UnknownSize', ...
                'Unknown thesis figure size preset: %s', char(sizeKey));
    end
end
