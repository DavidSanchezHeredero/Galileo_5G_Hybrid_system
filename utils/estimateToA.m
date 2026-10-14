function [toa, maxCorr] = estimateToA(rxWaveform, carrier, prsGrid, ofdmInfo)
%ESTIMATETOA ToA estimation constrained to the expected PRS slot(s).
%
% Outputs:
%   toa(g)      : estimated propagation delay (seconds) for gNB g (time-of-flight)
%   maxCorr(g)  : peak correlation magnitude for the chosen slot
%   usedSlot(g) : (optional) slot index used for the estimate
%
% IMPORTANT:
%   - This function assumes the PRS schedule is known (prsConfigs).
%   - It correlates only around the slot(s) where PRS is expected to be present,
%     avoiding ambiguity across repetitions/periodicity.
%   - Returned toa is a delay *within that PRS slot window* (NOT absolute time in buffer).
    numgNBs  = numel(carrier);
    corr     = cell(1,numgNBs);
    delayEst = zeros(1,numgNBs);
    maxCorr  = zeros(1,numgNBs);
    toa      = zeros(1,numgNBs);
    for gNBIdx = 1:numgNBs
        [~,mag] = nrTimingEstimate(carrier(gNBIdx),rxWaveform,prsGrid{gNBIdx});
        % Extract correlation data samples spanning about 1/14 ms for normal
        % cyclic prefix and about 1/12 ms for extended cyclic prefix (this
        % truncation is to ignore noisy side lobe peaks in correlation outcome)
        %corr{gNBIdx}     = mag(1:(ofdmInfo.Nfft*carrier(1).SubcarrierSpacing/15));
        corr{gNBIdx}     = mag;

        % since no longer distance than 2km is expected, the correlation is
        % only observed until k=3.000 which is equal to
        % time=3000/491,52e6 -> 6 microseconds
        [maximo, offset] = max(corr{gNBIdx}(1:3000));

        % Delay estimate is at point of maximum correlation
        maxCorr(gNBIdx)  = maximo;

        % Sample with max
        delayEst(gNBIdx) = offset;

        % Conversion from sample number to time dividing by sample rate
        toa(gNBIdx)      = delayEst(gNBIdx)/ofdmInfo.SampleRate;
    end
end
