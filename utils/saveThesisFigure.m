function saveThesisFigure(fig, fileName, fontAxis, fontTitle, fontLegend, gridColor, sizeSpec)
%SAVETHESISFIGURE Style and export a figure with thesis-ready dimensions.

    if nargin < 7
        sizeSpec = 'single';
    end

    outDir = fileparts(fileName);
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    styleThesisFigure(fig, fontAxis, fontTitle, fontLegend, gridColor, sizeSpec);

    [~, ~, ext] = fileparts(fileName);
    ext = lower(ext);
    switch ext
        case {'.pdf', '.svg', '.eps'}
            exportgraphics(fig, fileName, ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'white');
        case {'.png', '.jpg', '.jpeg', '.tif', '.tiff'}
            exportgraphics(fig, fileName, ...
                'Resolution', 400, ...
                'BackgroundColor', 'white');
        otherwise
            exportgraphics(fig, fileName, ...
                'Resolution', 400, ...
                'BackgroundColor', 'white');
    end
end
