function LoSplot(map3D,rxENU,trans,visibleCount,los)
    fontAxis   = 15;
    fontTitle  = 17;
    fontLegend = 13;
    gridColor  = [0.82 0.82 0.82];

    fig = createThesisFigure('three-d');
    show(map3D); hold on;
    
    % Receiver
    scatter3(rxENU(1), rxENU(2), rxENU(3), 100, 'k', 'filled');
    text(rxENU(1), rxENU(2), rxENU(3)+5, 'Receiver', 'Color', 'k', 'FontWeight', 'bold');
    
    % Satellites
    scatter3(trans(:,1), trans(:,2), trans(:,3), 80, 'g', 'filled');
    
    % Draw LoS lines.
    for s = 1:visibleCount
        tx(1) = trans(s,1);
        tx(2) = trans(s,2);
        tx(3) = trans(s,3);
        ptEnd = rxENU;
    
        % Use green for LoS and red for blocked links.
        lineColor = 'r';
        if  los(s)
            lineColor = 'g';
        end
        plot3([tx(1) ptEnd(1)], [tx(2) ptEnd(2)], [tx(3) ptEnd(3)], '-', ...
              'Color', lineColor, 'LineWidth', 2);
    end
    
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Line of sight between all sv and receiver.');
    axis equal; grid on; view(3);
    styleThesisFigure(fig, fontAxis, fontTitle, fontLegend, gridColor, 'three-d');
end
