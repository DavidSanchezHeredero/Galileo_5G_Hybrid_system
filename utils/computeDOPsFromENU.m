function dop = computeDOPsFromENU(txPositions_ENU, trueENU, visibilityMap, trueRange)
% computeDOPsFromENU
% Calculates GDOP, PDOP, HDOP, VDOP, TDOP per epoch from ENU positions.
%
% Input:
%   txPositions_ENU : [Ntx x 3 x Nep]  -> ENU svs positions
%   trueENU         : [Nep x 3]         -> ENU receiver positions per epoch
%   visibilityMap:  : [Ntx x Nep]
%
% Output:
%   dop: struct with fields (vectors [Nep x 1]):
%        .GDOP, .PDOP, .HDOP, .VDOP, .TDOP, .Ntransmitters
%
% Notes:
%   - Matrix G is built with rows [-uE, -uN, -uU, 1], where u is the
%     LoS unit vector from the user to the transmitter.
%   - If <4 satellites, NaN is returned.
    nr_case = false;
    [Ntx, dim, Nep] = size(txPositions_ENU);

    if Nep > 1
        assert(dim == 3, 'Second dimension of svs must be 3 (E,N,U).');
        assert(size(trueENU,1) == Nep && size(trueENU,2) == 3 && size(visibilityMap,2) == Nep , ...
            'trueENU must be [Nep x 3], and Nep must also match the number of epochs in txPositions_ENU and visibilityMap.');
        
        % Outputs
    GDOP = nan(Nep,1);
    PDOP = nan(Nep,1);
    HDOP = nan(Nep,1);
    VDOP = nan(Nep,1);
    TDOP = nan(Nep,1);
    Ntransmitters = zeros(Nep,1);

    % Inversion of the normal matrix.
    invSym = @(A) (A \ eye(size(A)));  % better performance than inv(A)

    for k = 1:Nep
        % Svs in epoch k
        tx_idx = find(visibilityMap(:,k) ==1);
        nlos    = numel(tx_idx);
        txk = txPositions_ENU(:,:,k);        % [Nsat x 3]
        usr     = trueENU(k,:);                      % [1 x 3]

        % LoS vector (user -> sat) ENU
        los            = zeros(Ntx,3)*NaN;
        los(tx_idx,:)  = txk(tx_idx,:) - usr;

        % Ranges and unit vector
        rho     = trueRange(:,k);               % [Nsat x 1]

        if isempty(rho)
            Ntransmitters(k) = 0;
            continue;
        end

        u = los(tx_idx,:) ./ rho(tx_idx);                          % unit vector (user -> sat)

        % Elevation filtering
        if nlos < 4
            % Four unknowns (x,y,z,dt) require at least four transmitters.
            continue;
        end

        % G matrix: rows [-uE, -uN, -uU, 1].
        G = [-u, ones(nlos,1)];

        % Q = (G'G)^{-1}
        % (if G'G is bad conditioned, \ handles the inversion better)
        GTG = G.' * G;
        % Conditioning check
        if rcond(GTG) < 1e-12
            % pseudo-inverse
            Q = pinv(GTG);
        else
            Q = invSym(GTG);
        end

        % DOPs
        qxx = Q(1,1); qyy = Q(2,2); qzz = Q(3,3); qtt = Q(4,4);
        GDOP(k) = sqrt(qxx + qyy + qzz + qtt);
        PDOP(k) = sqrt(qxx + qyy + qzz);
        HDOP(k) = sqrt(qxx + qyy);
        VDOP(k) = sqrt(qzz);
        TDOP(k) = sqrt(qtt);
    end
    elseif Nep == 1
        Nep = size(trueENU,1);
        Ntx = size(visibilityMap,2);

        % Outputs
        GDOP = nan(Nep,1);
        PDOP = nan(Nep,1);
        HDOP = nan(Nep,1);
        VDOP = nan(Nep,1);
        TDOP = nan(Nep,1);
        Ntransmitters = zeros(Nep,1);
    
        % Inversion of the normal matrix.
        invSym = @(A) (A \ eye(size(A)));  % better performance than inv(A)
    
        for k = 1:Nep
            % Svs in epoch k
            tx_idx = find(visibilityMap(k,:) ==1 );
            nlos    = numel(tx_idx);
            txk = txPositions_ENU;               % [Nbs x 3]

            usr     = trueENU(k,:);                      % [1 x 3]
    
            % LoS vector (user -> sat) ENU
            los            = zeros(Ntx,3)*NaN;
            los(tx_idx,:)  = txk(tx_idx,:) - usr;
    
            % Ranges and unit vector
            rho     = trueRange(:,k);               % [Nsat x 1]
    
            if isempty(rho)
                Ntransmitters(k) = 0;
                continue;
            end
    
            u = los(tx_idx,:) ./ rho(tx_idx);                          % unit vector (user -> sat)
    
            % Elevation filtering
            if nlos < 4
                % Four unknowns (x,y,z,dt) require at least four transmitters.
                continue;
            end
    
            % G matrix: rows [-uE, -uN, -uU, 1].
            G = [-u, ones(nlos,1)];
    
            % Q = (G'G)^{-1}
            % (if G'G is bad conditioned, \ handles the inversion better)
            GTG = G.' * G;
            % Conditioning check
            if rcond(GTG) < 1e-12
                % pseudo-inverse
                Q = pinv(GTG);
            else
                Q = invSym(GTG);
            end
    
            % DOPs
            qxx = Q(1,1); qyy = Q(2,2); qzz = Q(3,3); qtt = Q(4,4);
            GDOP(k) = sqrt(qxx + qyy + qzz + qtt);
            PDOP(k) = sqrt(qxx + qyy + qzz);
            HDOP(k) = sqrt(qxx + qyy);
            VDOP(k) = sqrt(qzz);
            TDOP(k) = sqrt(qtt);
        end
    end
    

    dop = struct('GDOP',GDOP,'PDOP',PDOP,'HDOP',HDOP,'VDOP',VDOP,'TDOP',TDOP,'Ntransmitters',Ntransmitters);
end
