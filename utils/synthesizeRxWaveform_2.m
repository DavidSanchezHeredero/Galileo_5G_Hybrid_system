function rxWaveform = synthesizeRxWaveform_2(txWaveforms, gNBPos, UEPos, ofdmInfo, fc, plCfg, addNoise, SNRdB)
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
    rxWaveform = zeros(size(txWaveforms,1) + maxDelay, 1);

    for g = 1:numgNBs
        % Path loss (LOS assumed here; if you want NLOS, pass flag)
        losFlag = true;
        PLdB = nrPathLoss(plCfg, fc, losFlag, gNBPos(g,:)', UEPos(:));
        PL = 10^(PLdB/10);
        
        %%%%
        symblen = ofdmInfo.SymbolLengths(2);
        cplen=ofdmInfo.CyclicPrefixLengths(2);
        cpdiff = ofdmInfo.CyclicPrefixLengths(1)-cplen;
        lengthHalfsubframe = sum(ofdmInfo.SymbolLengths)/2;
        deWav1=reshape(txWaveforms(:,g),lengthHalfsubframe,10);
        deWav2=deWav1(cpdiff+1:245760,:);
        deWav3=reshape(deWav2,symblen,1120/2);
        deWav4=fft(deWav3(cplen+1:symblen,:),4096,1);
        delShift = exp(-1i*2*pi*(1:ofdmInfo.Nfft)'*120e3*dFrac(g)/ofdmInfo.SampleRate);
        shiftWav = deWav4.*delShift;
        txWav1 = ifft(shiftWav);
        txWav2 = txWav1([4096-cplen+1:4096 1:4096],:);
        txWav3 = reshape(txWav2(:),lengthHalfsubframe-cpdiff,[]);
        txWav4 = txWav3([symblen-cplen-cpdiff+1:symblen-cplen 1:end],:);
        sig2 = [zeros(dInt(g),1); txWav4(:)];
        % Pad to same length and attenuate
        sig = [sig2; zeros(length(rxWaveform) - length(sig2), 1)] / sqrt(PL);

        rxWaveform = rxWaveform + sig;
    end


    % if addNoise
    %     rxWaveform = awgn(rxWaveform, SNRdB, 'measured');
    % end
end
