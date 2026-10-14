function rxWaveform = synthesizeRxWaveform(txWaveform, gNBPos, UEPos, ofdmInfo, fc, plCfg, addNoise, SNRdB)
    c = physconst('LightSpeed');
    numgNBs = size(gNBPos,1);

    % Geometric ToA
    diffPos  = gNBPos - UEPos;
    distance = vecnorm(diffPos,2,2);
    rtoa     = distance / c;

    % Fractional sample delays
    dSamp = rtoa * ofdmInfo.SampleRate;     % may be fractional
    dInt  = floor(dSamp);
    dFrac = dSamp - dInt;

    maxDelay = max(dInt) + 2;               % + margin
    rxWaveform = zeros(length(txWaveform{1}) + maxDelay, 1);

    % Fractional delay object (DSP toolbox)
    vfd = dsp.VariableFractionalDelay;

    for g = 1:numgNBs
        % Path loss (LOS assumed here; if you want NLOS, pass flag)
        losFlag = true;
        PLdB = nrPathLoss(plCfg, fc, losFlag, gNBPos(g,:)', UEPos(:));
        PL = 10^(PLdB/10);

        % Apply integer delay (pad)
        sig = [zeros(dInt(g),1); txWaveform{g}];

        % Apply fractional delay
        sig = vfd(sig, dFrac(g));

        % Pad to same length and attenuate
        sig = [sig; zeros(length(rxWaveform) - length(sig), 1)] / sqrt(PL);

        rxWaveform = rxWaveform + sig;
    end

    if addNoise
        rxWaveform = awgn(rxWaveform, SNRdB, 'measured');
    end
end
