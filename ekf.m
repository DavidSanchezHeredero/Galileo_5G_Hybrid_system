% =========================================================================
% Author: David Sanchez-Heredero Martinez
% Date: 2025-04-28
% Project: EKF Simulation
% =========================================================================
clc; clear; close all;

addpath('./utils','./config');
run('config');
load(rxENUFile, 'receiverPositions', 'samplingRate');
load(galLoSFile,'measPseudorange','Fs_gnss','svPositions_ENU');
load(nrObsFile,'tdoaEstim','tdoaError');
load(nrSimFile, 'Ts5g_pos');
load(nrLoSFile,'LoSbStationsENU','LoSMatrix');

epochs = size(receiverPositions,1);
dt      = 1/samplingRate;
dtGNSS  = 1/Fs_gnss;
dt5g    = Ts5g_pos;
losmx=LoSMatrix(1:dt5g/dt:end,:);

% Process noise
Q = blkdiag(1e-2*eye(3),1e-1*eye(3),10);
% Measurement noise variances
sigGNSS = 10;  Rgnss0 = sigGNSS^2;
sig5g   =  5;  R5g0    = sig5g^2;

% State transition (constant velocity + bias)
Fdt = @(dt) [eye(3), dt*eye(3), zeros(3,1);
             zeros(3), eye(3), zeros(3,1);
             zeros(1,6), 1];

x = zeros(7,1); x(1:3)=receiverPositions(1,1:3)';
P = diag([1e6*ones(3,1);100*ones(3,1);1e4]);
x_est = nan(7,epochs);

gnssCount = 0; nrCount = 0;

for epoch = 1:epochs
    t = dt*epoch;
    % 1. Predict step with correct dt
    if epoch>1
        dt_step = dt;
        F = Fdt(dt_step);
        x = F*x;
        P = F*P*F' + Q;
    end

    % 2. GNSS update
    if mod(t,dtGNSS) == 0 % true when dt is multiple of dt_GNSS
        gnssEpoch = t / dtGNSS;
        satIdx = find(~isnan(measPseudorange(:,(gnssEpoch))));
        numSv  = numel(satIdx);
        if numSv>0
            z = measPseudorange(satIdx,gnssEpoch);
            svPos = svPositions_ENU(satIdx,:,gnssEpoch);
            H = zeros(numSv,7); h = zeros(numSv,1);
            for i = 1:numSv
                ecd = norm(x(1:3)-svPos(i,:)');
                H(i,1:3)= (x(1:3)-svPos(i,:)')/ecd;
                H(i,7)=1;
                h(i) = ecd + x(7);
            end
            y = z - h;
            R = Rgnss0 * eye(numSv);

            % validation check
            if any(~isfinite(z)) || any(~isfinite(h)) || any(~isfinite(H(:))) || any(~isfinite(P(:)))
                disp("---- DEBUG GNSS ----");
                disp(["epoch=", epoch, " gnssEpoch=", gnssEpoch, " numSv=", numSv]);
                disp(["anyNaN svPos? ", any(~isfinite(svPos(:)))]);
                disp(["anyNaN H? ", any(~isfinite(H(:)))]);
                disp(["anyNaN h? ", any(~isfinite(h))]);
                disp(["anyNaN z? ", any(~isfinite(z))]);
                keyboard
            end

            S = H*P*H' + R;
            K = P*H'/S;
            x = x + K*y;
            P = (eye(7)-K*H)*P;
        end
    end

    % 3. 5G TDOA update
    if ~mod(t,dt5g) % true when dt is multiple of dt_5G
        nrEpoch = t / dt5g;
        vBS = find(losmx(nrEpoch,:));
        %vBS = find(~isnan(tdoaEstim(nrEpoch,:)));
        numBS = numel(vBS);
        if numBS >= 4
            z       = tdoaEstim(nrEpoch,vBS(2:end));
            bsPos   = LoSbStationsENU(vBS,:);
            numZ    = numBS-1;
            H = zeros(numZ,7); h = zeros(numZ,1);
            for i = 1:numZ
                d1 = norm(x(1:3) - bsPos(i,:)');
                d2 = norm(x(1:3) - bsPos(i+1,:)');
                h(i) = (d1 - d2)/physconst('LightSpeed');
                H(i,1:3) = (x(1:3)-bsPos(i,:)')'/d1 - (x(1:3)-bsPos(i+1,:)')'/d2;
            end
            y = z' - h;
            R = R5g0 * eye(numZ);
            % Validation check
            if any(~isfinite(z)) || any(~isfinite(h)) || any(~isfinite(H(:))) || any(~isfinite(P(:)))
                disp("---- DEBUG 5G ----");
                disp(["epoch=", epoch, " nrEpoch=", nrEpoch, " numBS=", numBS]);
                disp(["anyNaN bsPos? ", any(~isfinite(bsPos(:)))]);
                disp(["anyNaN H? ", any(~isfinite(H(:)))]);
                disp(["anyNaN h? ", any(~isfinite(h))]);
                disp(["anyNaN z? ", any(~isfinite(z))]);
                keyboard
            end

            S = H*P*H' + R;
            K = P*H'/S;
            x = x + K*y;
            P = (eye(7)-K*H)*P;
        end
    end

    x_est(:,epoch) = x;
end

ekf_error = receiverPositions(:,:) - x_est(1:3,:)';
for i = 1:size(receiverPositions,1)
    ekf_error(i) = norm(receiverPositions(i,:) - x_est(1:3,i)');
end
save(ekfSimFile,"x_est","ekf_error");
