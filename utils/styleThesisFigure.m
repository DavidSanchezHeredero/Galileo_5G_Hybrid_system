function styleThesisFigure(figs, fontAxis, fontTitle, fontLegend, gridColor, sizeSpec)
%STYLETHESISFIGURE Apply a consistent thesis-ready style to figures.

    if nargin < 6
        sizeSpec = 'single';
    end

    if ~iscell(sizeSpec)
        sizeSpec = repmat({sizeSpec}, 1, numel(figs));
    end

    for idxFig = 1:numel(figs)
        fig = figs(idxFig);
        figSizeCm = thesisFigureSize(sizeSpec{min(idxFig, numel(sizeSpec))});

        set(fig, ...
            'Color', 'w', ...
            'Units', 'centimeters', ...
            'Position', [2 2 figSizeCm], ...
            'PaperUnits', 'centimeters', ...
            'PaperPosition', [0 0 figSizeCm], ...
            'PaperSize', figSizeCm, ...
            'PaperPositionMode', 'auto', ...
            'InvertHardcopy', 'off');

        axList = [findall(fig, 'Type', 'axes'); findall(fig, 'Type', 'geoaxes')];
        axList = axList(isgraphics(axList));
        for k = 1:numel(axList)
            ax = axList(k);

            if isprop(ax, 'FontSize')
                ax.FontSize = fontAxis;
            end
            if isprop(ax, 'LineWidth')
                ax.LineWidth = 1.0;
            end
            if isprop(ax, 'Color')
                ax.Color = 'w';
            end
            if isprop(ax, 'GridColor')
                ax.GridColor = gridColor;
            end
            if isprop(ax, 'GridAlpha')
                ax.GridAlpha = 0.8;
            end
            if isprop(ax, 'XColor')
                ax.XColor = 'k';
            end
            if isprop(ax, 'YColor')
                ax.YColor = 'k';
            end
            if isprop(ax, 'ZColor')
                ax.ZColor = 'k';
            end
            if isprop(ax, 'FontWeight')
                ax.FontWeight = 'bold';
            end
            if isprop(ax, 'LatitudeAxis') && isprop(ax.LatitudeAxis, 'Color')
                ax.LatitudeAxis.Color = 'k';
            end
            if isprop(ax, 'LongitudeAxis') && isprop(ax.LongitudeAxis, 'Color')
                ax.LongitudeAxis.Color = 'k';
            end

            try
                box(ax, 'on');
            catch
            end

            styleAxisText(ax, 'XLabel', fontAxis);
            styleAxisText(ax, 'YLabel', fontAxis);
            styleAxisText(ax, 'ZLabel', fontAxis);
            styleAxisText(ax, 'Title', fontTitle);
        end

        lgdList = findall(fig, 'Type', 'legend');
        for k = 1:numel(lgdList)
            lgd = lgdList(k);
            lgd.FontSize = fontLegend;
            if isprop(lgd, 'Color')
                lgd.Color = 'w';
            end
            if isprop(lgd, 'TextColor')
                lgd.TextColor = 'k';
            end
        end
    end
end

function styleAxisText(ax, propertyName, fontSize)
    if ~isprop(ax, propertyName)
        return;
    end

    txt = ax.(propertyName);
    if ~isgraphics(txt)
        return;
    end

    if isprop(txt, 'FontSize')
        txt.FontSize = fontSize;
    end
    if isprop(txt, 'FontWeight')
        txt.FontWeight = 'bold';
    end
    if isprop(txt, 'Color')
        txt.Color = 'k';
    end
end
