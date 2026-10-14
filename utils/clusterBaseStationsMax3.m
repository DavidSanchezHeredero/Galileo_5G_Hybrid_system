function groups = clusterBaseStationsMax3(LoSbStationsENU, sGroup)
    numBS       = length(LoSbStationsENU(:,1));
    assigned    = false(numBS, 1);
    groups      = {};
    groupIdx    = 1;

    % Copy XY coordinates; Z is not used for clustering.
    positions = LoSbStationsENU(:, 1:2);

    % Continue until every base station has been assigned.
    while any(~assigned)
        idx = find(~assigned, 1);
        center = positions(idx, :);

        % Compute distances from the current base station to all unassigned stations.
        availableIdx = find(~assigned);
        availableCoords = positions(availableIdx, :);

        % Sort by distance to the current point.
        dists = vecnorm(availableCoords - center, 2, 2);
        [~, sortedRelIdx] = sort(dists);

        % Select up to sGroup nearest stations.
        selectedRel = sortedRelIdx(1:min(sGroup, length(sortedRelIdx)));
        selectedIdx = availableIdx(selectedRel);

        % Assign the selected stations to the current group.
        groups{groupIdx} = selectedIdx;
        assigned(selectedIdx) = true;
        groupIdx = groupIdx + 1;
    end
end
