function rxWaveform = synthesizeRxWaveform_2(txWaveforms, gNBPos, UEPos, ofdmInfo, fc, plCfg, addNoise, SNRdB)
    c = physconst('LightSpeed');
    numgNBs = size(gNBPos,1);
    txSamples = size(txWaveforms,1);

    % Geometric ToA
    diffPos  = gNBPos - UEPos;
    distance = vecnorm(diffPos,2,2);
    rtoa     = distance / c;

    % Fractional sample delays
    dSamp = rtoa * ofdmInfo.SampleRate;     % may be fractional
    dInt  = floor(dSamp);
    dFrac = gpuArray(single(dSamp - dInt));

    maxDelay = max(dInt) + 2;               % + margin
    rxWaveform = zeros(size(txWaveforms,1) + maxDelay, 1);
    
    % t created directly as gpuArray to force later operations in GPU
    t = gpuArray(single(1:txSamples))'/txSamples;
 
    for g = 1:numgNBs
        % Path loss (LOS assumed here; if you want NLOS, pass flag)
        losFlag = true;
        PLdB = nrPathLoss(plCfg, fc, losFlag, gNBPos(g,:)', UEPos(:));
        PL = 10^(PLdB/10);

        % Convert full signal to frequency domain
        dftWav = fft(txWaveforms(:,g));

        % Create fractional delay shift matrix
        dftShift = exp(-1i*2*pi*t*dFrac(g));
        
        % Apply fractional delay shift
        dshiftWav = dftWav.*dftShift;

        % Return to time domain
        fractDelayedWave = ifft(dshiftWav);

        % Apply integer delay   
        delayedWave = [zeros(dInt(g),1); fractDelayedWave];

        % Pad to same length and attenuate
        delayedWave = [delayedWave; zeros(length(rxWaveform) - length(delayedWave), 1)] / sqrt(PL);
        
        % Add to total received wave
        rxWaveform = rxWaveform + delayedWave;
    end

    % if addNoise
    %     rxWaveform = awgn(rxWaveform, SNRdB, 'measured');
    % end
end
