function report = Hybrid_EKF_Publication_Export(results, settings, outputDir)
% HYBRID_EKF_PUBLICATION_EXPORT  Reproducible Springer Nature outputs.
%   report = Hybrid_EKF_Publication_Export(results,settings,outputDir)
%   Inputs: results/settings produced by Hybrid_EKF_Integrity.m(.mlx).
%   To regenerate WITHOUT rerunning the EKF:
%     S = load('Hybrid_EKF_Integrity_results.mat','results','settings');
%     report = Hybrid_EKF_Publication_Export(S.results,S.settings,'springer_export');
%
%   Outputs: CSV (per-epoch and summary), LaTeX booktabs tables, vector PDF,
%   600 dpi PNG, optional editable MATLAB FIG, and a provenance manifest.
%   No standalone GNSS/5G metrics are inferred from the hybrid estimator.
%   All statistics are descriptive of ONE simulation (time-correlated epochs).
%   Neither a confidence-level coverage claim nor certified PL is implied.
%
%   Springer's official journal-specific author instructions take precedence.

    if nargin < 3 || isempty(outputDir)
        outputDir = fullfile(pwd,'springer_export');
    end
    if isstring(outputDir), outputDir = char(outputDir); end
    assert(isstruct(results) && isstruct(settings), ...
        'Expected results and settings structs from Hybrid_EKF_Integrity.');
    assert(exist('exportgraphics','file') ~= 0, ...
        'This exporter needs MATLAB R2020a or later (exportgraphics).');

    figDir = fullfile(outputDir,'figures');
    tabDir = fullfile(outputDir,'tables');
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    if ~exist(figDir,'dir'), mkdir(figDir); end
    if ~exist(tabDir,'dir'), mkdir(tabDir); end

    t = asRow(results.time_s);
    N = numel(t);
    required = {'errorENU','horizontalError','verticalError', ...
                'HPL0','VPL0','HPL','VPL','NEESpos', ...
                'GNSSaccepted','NRaccepted','GNSSalarm','NRalarm', ...
                'NISgnss','NIS5g','NISThreshGnss','NISThresh5g', ...
                'numberGNSS','numberTDOA','Ppos'};
    for i = 1:numel(required)
        assert(isfield(results,required{i}), ...
            'Missing results.%s. Rerun the up-to-date EKF.', required{i});
    end
    assert(isequal(size(results.Ppos),[3 3 N]), 'Expected Ppos [3 x 3 x N].');
    assert(isequal(size(results.errorENU),[3 N]), 'Expected errorENU [3 x N].');
    for i = 1:numel(required)
        v = results.(required{i});
        if isvector(v), assert(numel(v)==N,'Wrong length of results.%s.',required{i}); end
    end
    assert(isfield(settings,'burnInEpochs') && isfield(settings,'HAL_m') && ...
           isfield(settings,'VAL_m'), 'Missing run configuration.');
    assert(N > settings.burnInEpochs && settings.burnInEpochs>=0);

    % The evaluation subset is explicit; prediction-only epochs are retained.
    assessed = false(1,N);
    assessed(settings.burnInEpochs+1:end) = true;
    eh = asRow(results.horizontalError);
    ev = asRow(results.verticalError);
    hp = asRow(results.HPL);
    vp = asRow(results.VPL);
    hp0 = asRow(results.HPL0);
    vp0 = asRow(results.VPL0);
    nG = asRow(results.numberGNSS);
    n5 = asRow(results.numberTDOA);
    acceptedG = logical(asRow(results.GNSSaccepted));
    accepted5 = logical(asRow(results.NRaccepted));
    alarmG = logical(asRow(results.GNSSalarm));
    alarm5 = logical(asRow(results.NRalarm));
    alarmAny = alarmG | alarm5;
    nisG = asRow(results.NISgnss);
    nis5 = asRow(results.NIS5g);
    thrG = asRow(results.NISThreshGnss);
    thr5 = asRow(results.NISThresh5g);
    nees = asRow(results.NEESpos);

    % P must be positive definite; otherwise a nominal Gaussian PL is invalid.
    pdCov = false(1,N);
    sigmaPos = nan(3,N);
    for k=1:N
        pk=(results.Ppos(:,:,k)+results.Ppos(:,:,k)')/2;
        [~,flag]=chol(pk);
        pdCov(k)=flag==0 && all(isfinite(pk(:)));
        if pdCov(k), sigmaPos(:,k)=sqrt(max(diag(pk),0)); end
    end
    valid = assessed & pdCov & all(isfinite(results.errorENU),1) & ...
        isfinite(eh) & isfinite(ev) & isfinite(hp) & isfinite(vp) & ...
        isfinite(hp0) & isfinite(vp0) & hp>0 & vp>0;
    nAssessed=nnz(assessed);
    nValid=nnz(valid);
    assert(nValid>0, 'No finite, positive-definite covariance/PL epochs after burn-in.');

    noAlert = valid & ~alarmAny;
    admitted = noAlert & hp<=settings.HAL_m & vp<=settings.VAL_m;
    bothAccepted = valid & acceptedG & accepted5;
    onlyG = valid & acceptedG & ~accepted5;
    only5 = valid & ~acceptedG & accepted5;
    predictionOnly = valid & ~acceptedG & ~accepted5;
    exceedH = valid & eh>hp;
    exceedV = valid & ev>vp;
    hazardH = admitted & eh>settings.HAL_m;
    hazardV = admitted & ev>settings.VAL_m;
    hazardAny = hazardH | hazardV;
    testedG = assessed & isfinite(nisG) & isfinite(thrG) & thrG>0;
    tested5 = assessed & isfinite(nis5) & isfinite(thr5) & thr5>0;
    testedNEES = valid & isfinite(nees);

    if isfield(results,'summary') && isfield(results.summary,'referenceIndexSupplied')
        referenceProvided = logical(results.summary.referenceIndexSupplied);
    else
        referenceProvided = false;
    end
    if isfield(settings,'initialPositionFromTruth')
        truthInitialization = logical(settings.initialPositionFromTruth);
    else
        truthInitialization = true; % Treat missing provenance conservatively.
    end
    if isfield(settings,'gateMeasurements')
        gate = logical(settings.gateMeasurements);
    else
        gate = false;
    end
    if isfield(results,'fault') && isfield(results.fault,'enabled')
        faultEnabled = logical(results.fault.enabled);
    else
        faultEnabled = false;
    end
    if isfield(settings,'maxGNSSBias_m') && isfield(settings,'maxTDOABias_m')
        boundedMode = settings.maxGNSSBias_m>0 || settings.maxTDOABias_m>0;
    else
        boundedMode = false;
    end
    if ~referenceProvided
        status='PROVISIONAL: TDOA reference not independently verified';
    elseif truthInitialization
        status='RESEARCH ONLY: filter initialised using simulation truth';
    else
        status='RESEARCH ONLY: model-conditional PL, not certified integrity';
    end

    % DESCRIPTIVE results. No independence assumption is used in these ratios.
    m = struct();
    m.nAssessed=nAssessed; m.nValid=nValid;
    m.nPredictionOnly=nnz(predictionOnly);
    m.nOnlyGNSS=nnz(onlyG); m.nOnly5G=nnz(only5);
    m.nBoth=nnz(bothAccepted);
    m.nAnyUpdate=nnz(valid & (acceptedG | accepted5));
    m.nGNSSUpdates=nnz(assessed & acceptedG);
    m.n5GUpdates=nnz(assessed & accepted5);
    m.nGNSSTests=nnz(testedG); m.n5GTests=nnz(tested5);
    m.nGNSSAlarms=nnz(assessed & alarmG);
    m.n5GAlarms=nnz(assessed & alarm5);
    m.nAlarmEpochs=nnz(valid & alarmAny);
    m.nAdmitted=nnz(admitted);
    m.computability_percent=100*nValid/nAssessed;
    m.thresholdAdmission_percent=100*m.nAdmitted/nAssessed;
    m.horizRMSE_m=sqrt(mean(eh(valid).^2));
    m.vertRMSE_m=sqrt(mean(ev(valid).^2));
    m.horizP95Error_m=localPercentile(eh(valid),95);
    m.vertP95Error_m=localPercentile(ev(valid),95);
    m.horizMedianPL_m=localPercentile(hp(valid),50);
    m.vertMedianPL_m=localPercentile(vp(valid),50);
    m.horizP95PL_m=localPercentile(hp(valid),95);
    m.vertP95PL_m=localPercentile(vp(valid),95);
    m.horizExceedances=nnz(exceedH);
    m.vertExceedances=nnz(exceedV);
    m.horizContainment_percent=100*(1-m.horizExceedances/nValid);
    m.vertContainment_percent=100*(1-m.vertExceedances/nValid);
    m.horizHazardCount=nnz(hazardH);
    m.vertHazardCount=nnz(hazardV);
    m.anyHazardCount=nnz(hazardAny);
    m.gnssNISAlarm_percent=localRate(nisG(testedG)>thrG(testedG));
    m.fiveGNISAlarm_percent=localRate(nis5(tested5)>thr5(tested5));
    m.meanNEES=mean(nees(testedNEES));
    m.nNEES=nnz(testedNEES);
    chiLower=2*gammaincinv(0.025,3/2);
    chiUpper=2*gammaincinv(0.975,3/2);
    m.neesOutside95_percent=localRate(nees(testedNEES)<chiLower | ...
                                          nees(testedNEES)>chiUpper);
    m.neesUpper95=chiUpper;
    m.neesLower95=chiLower;

    % Per-epoch machine-readable results allow auditing every table entry.
    state = zeros(1,N);
    state(acceptedG & ~accepted5)=1;
    state(~acceptedG & accepted5)=2;
    state(acceptedG & accepted5)=3;
    state(alarmAny)=4; % Alarm state takes priority over update classification.
    epochTable = table((1:N)',t',assessed',valid', ...
        results.errorENU(1,:)',results.errorENU(2,:)',results.errorENU(3,:)', ...
        eh',ev',hp0',vp0',hp',vp', ...
        sigmaPos(1,:)',sigmaPos(2,:)',sigmaPos(3,:)', ...
        nG',n5',acceptedG',accepted5',alarmG',alarm5', ...
        nisG',thrG',nis5',thr5',nees', ...
        admitted',exceedH',exceedV',hazardAny',state', ...
        'VariableNames',{'epoch','time_s','assessed','valid', ...
        'errorE_m','errorN_m','errorU_m','errorH_m','errorV_m', ...
        'HPL0_m','VPL0_m','HPL_m','VPL_m', ...
        'stdE_m','stdN_m','stdU_m','numGNSS','numTDOA', ...
        'GNSSaccepted','NRaccepted','GNSSalarm','NRalarm', ...
        'NISgnss','thresholdGNSS','NIS5g','threshold5g','NEESpos', ...
        'admittedIllustrative','exceedHPL','exceedVPL', ...
        'observedHazardNoAlert','updateCategory'});
    writetable(epochTable,fullfile(outputDir,'hybrid_ekf_epoch_results.csv'));
    names=fieldnames(m);
    values=cellfun(@(f) m.(f),names);
    metricsTable=table(names,values,'VariableNames',{'metric','value'});
    writetable(metricsTable,fullfile(outputDir,'hybrid_ekf_metrics.csv'));
    report=struct('metrics',m,'modelStatus',status,'referenceProvided', ...
        referenceProvided,'truthInitialization',truthInitialization, ...
        'faultEnabled',faultEnabled,'gatingEnabled',gate, ...
        'validMask',valid,'assessmentMask',assessed, ...
        'admittedMask',admitted,'outputDir',outputDir);
    save(fullfile(outputDir,'hybrid_ekf_publication_report.mat'), ...
         'report','settings','-v7');

    % LaTeX-ready booktabs tables. Figure titles/captions stay in the .tex.
    writeAvailabilityTable(fullfile(tabDir,'table_ekf_availability.tex'),m,settings);
    writeAccuracyTable(fullfile(tabDir,'table_ekf_accuracy.tex'),m,settings,boundedMode);
    writeConsistencyTable(fullfile(tabDir,'table_ekf_consistency.tex'),m,settings);
    % Stamp individual tables so that copy/paste cannot silently erase key caveats.
    stampTableStatus(fullfile(tabDir,'table_ekf_availability.tex'),referenceProvided,truthInitialization);
    stampTableStatus(fullfile(tabDir,'table_ekf_accuracy.tex'),referenceProvided,truthInitialization);
    stampTableStatus(fullfile(tabDir,'table_ekf_consistency.tex'),referenceProvided,truthInitialization);
    writeFigureSnippets(fullfile(outputDir,'figures_ekf.tex'));
    writeInclusionSnippets(fullfile(outputDir,'insert_into_manuscript.tex'));
    writeManifest(fullfile(outputDir,'README_provenance.txt'),status,settings,m, ...
                  referenceProvided,truthInitialization,gate,faultEnabled,boundedMode,results);

    % Publication artwork: standard 17.8-cm double-column width; captions
    % and figure numbers are defined in the manuscript, NOT inside the plots.
    mask = valid;
    xf = t(mask);
    c.Black=[0.10 0.10 0.10];
    c.Blue=[0.00 0.34 0.63];
    c.Teal=[0.00 0.48 0.45];
    c.Red=[0.73 0.20 0.14];
    c.Grey=[0.45 0.45 0.45];

    f=makeFigure(17.8,8.5);
    plot(xf,eh(mask),'Color',c.Black,'LineWidth',1.0,'DisplayName','Horizontal error'); hold on;
    plot(xf,hp0(mask),'--','Color',c.Blue,'LineWidth',1.1,'DisplayName','Nominal HPL');
    if boundedMode
        plot(xf,hp(mask),'-','Color',c.Teal,'LineWidth',1.2,'DisplayName','Bias-bounded HPL');
    end
    yline(settings.HAL_m,':','Color',c.Red,'LineWidth',1.0, ...
          'DisplayName','Illustrative HAL');
    xlabel('Time (s)'); ylabel('Horizontal error / bound (m)');
    formatAxes(f); exportFigure(f,figDir,'fig_ekf_horizontal');

    f=makeFigure(17.8,8.5);
    plot(xf,ev(mask),'Color',c.Black,'LineWidth',1.0,'DisplayName','Vertical error'); hold on;
    plot(xf,vp0(mask),'--','Color',c.Blue,'LineWidth',1.1,'DisplayName','Nominal VPL');
    if boundedMode
        plot(xf,vp(mask),'-','Color',c.Teal,'LineWidth',1.2,'DisplayName','Bias-bounded VPL');
    end
    yline(settings.VAL_m,':','Color',c.Red,'LineWidth',1.0, ...
          'DisplayName','Illustrative VAL');
    xlabel('Time (s)'); ylabel('Absolute vertical error / bound (m)');
    formatAxes(f); exportFigure(f,figDir,'fig_ekf_vertical');

    f=makeFigure(17.8,14.5);
    tl=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
    coords={'East','North','Up'};
    for i=1:3
        ax=nexttile(tl);
        e=results.errorENU(i,mask);
        ss=sigmaPos(i,mask);
        plot(ax,xf,e,'Color',c.Black,'LineWidth',0.95,'DisplayName','Estimation error'); hold(ax,'on');
        plot(ax,xf,ss,'--','Color',c.Blue,'LineWidth',0.95,'DisplayName','+1 SD');
        plot(ax,xf,-ss,'--','Color',c.Blue,'LineWidth',0.95, ...
             'HandleVisibility','off');
        ylabel(ax,sprintf('%s (m)',coords{i})); grid(ax,'on');
        if i<3, ax.XTickLabel=[]; end
        if i==1, legend(ax,'Location','best','Box','off'); end
    end
    xlabel(tl,'Time (s)');
    formatAxes(f,false); exportFigure(f,figDir,'fig_ekf_enu');

    f=makeFigure(17.8,10.4);
    tl=tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl);
    rG=nan(1,N); rG(testedG)=nisG(testedG)./thrG(testedG);
    plot(ax,t(assessed),rG(assessed),'o-','Color',c.Blue, ...
         'MarkerSize',2,'LineWidth',0.9); hold(ax,'on');
    yline(ax,1,'--','Color',c.Red,'LineWidth',1.0);
    ylabel(ax,'GNSS NIS / threshold'); grid(ax,'on'); ax.XTickLabel=[];
    ax=nexttile(tl);
    r5=nan(1,N); r5(tested5)=nis5(tested5)./thr5(tested5);
    plot(ax,t(assessed),r5(assessed),'o-','Color',c.Teal, ...
         'MarkerSize',2,'LineWidth',0.9); hold(ax,'on');
    yline(ax,1,'--','Color',c.Red,'LineWidth',1.0);
    ylabel(ax,'5G NIS / threshold'); grid(ax,'on');
    xlabel(tl,'Time (s)');
    formatAxes(f,false); exportFigure(f,figDir,'fig_ekf_nis');

    f=makeFigure(17.8,8.5);
    plot(t(testedNEES),nees(testedNEES),'Color',c.Black,'LineWidth',0.95); hold on;
    yline(chiLower,'--','Color',c.Blue,'LineWidth',1.1, ...
        'DisplayName','2.5% quantile');
    yline(chiUpper,'-.','Color',c.Red,'LineWidth',1.1, ...
        'DisplayName','97.5% quantile');
    xlabel('Time (s)'); ylabel('Position NEES');
    legend('NEES','2.5% chi-square limit','97.5% chi-square limit', ...
           'Location','best','Box','off');
    formatAxes(f,false); exportFigure(f,figDir,'fig_ekf_nees');

    f=makeFigure(17.8,7.3);
    stairs(xf,state(mask),'Color',c.Black,'LineWidth',1.05);
    xlabel('Time (s)'); ylabel('Update and alarm state');
    yticks(0:4); yticklabels({'Predict','GNSS','5G','Both','Alarm'});
    ylim([-0.4 4.4]); grid on;
    formatAxes(f,false); exportFigure(f,figDir,'fig_ekf_updates');

    f=makeFigure(17.8,8.4);
    hClear=plot(hp(noAlert),eh(noAlert),'.','Color',c.Blue,'MarkerSize',7); hold on;
    legendHandles=hClear;
    legendLabels={'No contemporaneous NIS alarm'};
    if any(valid & alarmAny)
        hAlarm=plot(hp(valid & alarmAny),eh(valid & alarmAny),'x', ...
            'Color',c.Red,'MarkerSize',4,'LineWidth',0.8);
        legendHandles=[legendHandles,hAlarm]; %#ok<AGROW>
        legendLabels{end+1}='Alarm epoch'; %#ok<AGROW>
    end
    maxVal=max([hp(valid),eh(valid),1]);
    hDiag=plot([0 maxVal],[0 maxVal],'--','Color',c.Black,'LineWidth',1.0);
    xlabel('Horizontal PL (m)'); ylabel('Horizontal error (m)');
    legend([legendHandles,hDiag],[legendLabels,{'Error = PL'}], ...
           'Location','best','Box','off');
    axis square; xlim([0 maxVal]); ylim([0 maxVal]);
    formatAxes(f,false); exportFigure(f,figDir,'fig_ekf_pl_scatter');

    fprintf('\n[EKF journal export] %s\n',outputDir);
    fprintf('  Evaluated: %d; valid: %d; illustrative threshold admission: %.1f%%.\n', ...
        m.nAssessed,m.nValid,m.thresholdAdmission_percent);
    fprintf('  7 figures x (vector PDF, 600 dpi PNG) + LaTeX tables + CSV.\n');
    if ~referenceProvided
        warning('PUBLICATION PROVISIONAL: 5G reference indices are not independently verified.');
    end
    if truthInitialization
        warning('PUBLICATION PROVISIONAL: initial EKF position uses ground truth.');
    end
end

function v=asRow(v)
    v=v(:)';
end
function y=localRate(logicalVector)
    if isempty(logicalVector),y=NaN;else,y=100*mean(logicalVector);end
end
function x=localPercentile(data,p)
    data=sort(data(isfinite(data)));
    if isempty(data),x=NaN;return;end
    if numel(data)==1,x=data(1);return;end
    pos=1+(numel(data)-1)*p/100;
    low=floor(pos); hi=ceil(pos);
    x=data(low)+(pos-low)*(data(hi)-data(low));
end
function s=fmt(x)
    if ~isfinite(x),s='--';else,s=sprintf('%.2f',x);end
end
function s=intfmt(x)
    if ~isfinite(x),s='--';else,s=sprintf('%d',round(x));end
end
function writeTextLines(file,lines)
    fid=fopen(file,'w');
    assert(fid>0,'Could not write %s',file);
    clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
    for k=1:numel(lines)
        fprintf(fid,'%s\n',lines{k});
    end
end
function writeAvailabilityTable(file,m,settings)
    L={
        '\begin{table}[t]'
        '\centering'
        '\caption{Hybrid EKF update, alarm and threshold-admission statistics for one trajectory. Values exclude burn-in epochs. Limits are illustrative, not application certification criteria.}'
        '\label{tab:ekf-export-availability}'
        '\small'
        '\begin{tabular}{@{}lr@{}}'
        '\toprule'
        'Indicator & Hybrid EKF \\'
        '\midrule'
        sprintf('Evaluated epochs & %s \\',intfmt(m.nAssessed))
        sprintf('Valid covariance / PL epochs & %s \\',intfmt(m.nValid))
        sprintf('Prediction-only epochs & %s \\',intfmt(m.nPredictionOnly))
        sprintf('GNSS-only accepted epochs & %s \\',intfmt(m.nOnlyGNSS))
        sprintf('5G-only accepted epochs & %s \\',intfmt(m.nOnly5G))
        sprintf('Both modalities accepted & %s \\',intfmt(m.nBoth))
        sprintf('GNSS / 5G test events & %s / %s \\',intfmt(m.nGNSSTests),intfmt(m.n5GTests))
        sprintf('GNSS / 5G alarms & %s / %s \\',intfmt(m.nGNSSAlarms),intfmt(m.n5GAlarms))
        sprintf('At least one alarm (epochs) & %s \\',intfmt(m.nAlarmEpochs))
        sprintf('Computability (\%%) & %s \\',fmt(m.computability_percent))
        sprintf('Below illustrative HAL/VAL, no alarm (\%%) & %s \\',fmt(m.thresholdAdmission_percent))
        '\bottomrule'
        '\end{tabular}'
        sprintf('\\par\\smallskip{\\footnotesize Evaluated: $N=%d$ after burn-in. Illustrative $HAL=%s$ m, $VAL=%s$ m. Predictions are retained in valid outputs; temporal freshness is not certified.}',m.nAssessed,fmt(settings.HAL_m),fmt(settings.VAL_m))
        '\end{table}'
    };
    writeTextLines(file,L);
end
function writeAccuracyTable(file,m,settings,boundedMode)
    if boundedMode
        boundDesc='Nominal Gaussian plus assumed deterministic componentwise bias bounds';
    else
        boundDesc='Fault-free nominal Gaussian covariance bounds (zero deterministic bias allowances)';
    end
    L={
        '\begin{table}[t]'
        '\centering'
        sprintf('\\caption{Hybrid EKF accuracy and protection-bound diagnostics, %d valid epochs. %s.}',m.nValid,boundDesc)
        '\label{tab:ekf-export-accuracy}'
        '\small'
        '\begin{tabular}{@{}lrr@{}}'
        '\toprule'
        'Indicator & Horizontal & Vertical \\'
        '\midrule'
        sprintf('RMSE (m) & %s & %s \\',fmt(m.horizRMSE_m),fmt(m.vertRMSE_m))
        sprintf('95th-percentile error (m) & %s & %s \\',fmt(m.horizP95Error_m),fmt(m.vertP95Error_m))
        sprintf('Median PL (m) & %s & %s \\',fmt(m.horizMedianPL_m),fmt(m.vertMedianPL_m))
        sprintf('95th-percentile PL (m) & %s & %s \\',fmt(m.horizP95PL_m),fmt(m.vertP95PL_m))
        sprintf('PL exceedances (count) & %s & %s \\',intfmt(m.horizExceedances),intfmt(m.vertExceedances))
        sprintf('Empirical containment (\%%) & %s & %s \\',fmt(m.horizContainment_percent),fmt(m.vertContainment_percent))
        sprintf('Observed hazardous no-alert events (count) & %s & %s \\',intfmt(m.horizHazardCount),intfmt(m.vertHazardCount))
        '\bottomrule'
        '\end{tabular}'
        sprintf('\\par\\smallskip{\\footnotesize Error and PL are computed from the same %d evaluable filter epochs; illustrative $HAL=%s$ m and $VAL=%s$ m. Empirical frequencies on temporally correlated epochs do not validate a rare-event integrity probability.}',m.nValid,fmt(settings.HAL_m),fmt(settings.VAL_m))
        '\end{table}'
    };
    writeTextLines(file,L);
end
function writeConsistencyTable(file,m,settings)
    L={
        '\begin{table}[t]'
        '\centering'
        '\caption{Hybrid EKF innovation screening and retrospective covariance-consistency diagnostics. Events and epochs are temporally dependent.}'
        '\label{tab:ekf-export-consistency}'
        '\small'
        '\begin{tabular}{@{}lr@{}}'
        '\toprule'
        'Indicator & Observed value \\'
        '\midrule'
        sprintf('GNSS / 5G NIS test events & %s / %s \\',intfmt(m.nGNSSTests),intfmt(m.n5GTests))
        sprintf('GNSS NIS threshold exceedance (\%%) & %s \\',fmt(m.gnssNISAlarm_percent))
        sprintf('5G NIS threshold exceedance (\%%) & %s \\',fmt(m.fiveGNISAlarm_percent))
        sprintf('Position NEES evaluated epochs & %s \\',intfmt(m.nNEES))
        sprintf('Mean position NEES & %s \\',fmt(m.meanNEES))
        sprintf('NEES outside nominal 95\%% interval (\%%) & %s \\',fmt(m.neesOutside95_percent))
        '\bottomrule'
        '\end{tabular}'
        sprintf('\\par\\smallskip{\\footnotesize NIS event-level $P_{FA}=%s$. Marginal position NEES compared with $\\chi^2_3$ 2.5\%% and 97.5\%% quantiles. Gating changes the selected error distribution.}',sprintf('%.1e',settings.pFalseAlarm))
        '\end{table}'
    };
    writeTextLines(file,L);
end
function stampTableStatus(file,referenceProvided,truthInitialization)
    warnings={};
    if ~referenceProvided
        warnings{end+1}='5G TDOA reference association unverified'; %#ok<AGROW>
    end
    if truthInitialization
        warnings{end+1}='EKF initialized using simulation truth'; %#ok<AGROW>
    end
    if isempty(warnings),return;end
    fid=fopen(file,'r'); assert(fid>0,'Cannot read %s',file);
    txt=fread(fid,'*char')'; fclose(fid);
    mark='\end{table}';
    warnText=['\par\smallskip{\footnotesize \emph{Provisional: ' strjoin(warnings,'; ') '.}}'];
    assert(contains(txt,mark),'No LaTeX table environment in %s.',file);
    txt=strrep(txt,mark,[warnText newline mark]);
    fid=fopen(file,'w'); assert(fid>0,'Cannot write %s',file);
    fprintf(fid,'%s',txt); fclose(fid);
end
function writeFigureSnippets(file)
    figures={
        'fig_ekf_horizontal','Time evolution of hybrid EKF horizontal error and the Gaussian model-conditional HPL. The HAL is illustrative and not a validated application requirement.','fig:ekf-export-horizontal';
        'fig_ekf_vertical','Time evolution of absolute vertical EKF error and conditional VPL. The VAL is illustrative.','fig:ekf-export-vertical';
        'fig_ekf_enu','East, North and Up estimation errors with symmetric one-standard-deviation bounds from posterior covariance (not protection levels).','fig:ekf-export-enu';
        'fig_ekf_nis','GNSS and 5G innovation statistics normalized by the corresponding per-event chi-square alarm thresholds. Values exceeding unity trigger an alarm.','fig:ekf-export-nis';
        'fig_ekf_nees','Retrospective position NEES with the nominal 95 percent two-sided marginal chi-square interval (three degrees of freedom); samples along one trajectory are dependent.','fig:ekf-export-nees';
        'fig_ekf_updates','Hybrid filter status by epoch. Prediction-only, accepted measurement modalities and NIS alarm events are distinguished.','fig:ekf-export-updates';
        'fig_ekf_pl_scatter','Horizontal error versus model-conditional HPL, with error-equals-bound reference and distinct alarm/no-alarm epochs.','fig:ekf-export-scatter'
    };
    L={'% Generated by Hybrid_EKF_Publication_Export.m.'; ...
       '% Insert these figure environments where relevant; do not reuse standalone figure captions.'};
    for k=1:size(figures,1)
        lines={
            ''
            '\begin{figure*}[t]'
            '\centering'
            sprintf('\\includegraphics[width=\\textwidth]{springer_export/figures/%s.pdf}',figures{k,1})
            sprintf('\\caption{%s}',figures{k,2})
            sprintf('\\label{%s}',figures{k,3})
            '\end{figure*}'
        };
        L=[L;lines]; %#ok<AGROW>
    end
    writeTextLines(file,L);
end
function writeInclusionSnippets(file)
    L={
        '% Place near the hybrid experimental-results subsection of sn-article-integrity.tex.'
        '% Requires existing preamble: \usepackage{graphicx,booktabs}.'
        '% Compile with the publication output directory alongside the main .tex.'
        '\input{springer_export/tables/table_ekf_availability.tex}'
        '\input{springer_export/tables/table_ekf_accuracy.tex}'
        '\input{springer_export/tables/table_ekf_consistency.tex}'
        '\input{springer_export/figures_ekf.tex}'
    };
    writeTextLines(file,L);
end
function writeManifest(file,status,settings,m,referenceProvided,truthInit,gate,faultEnabled,boundedMode,results)
    L={
        'SPRINGER NATURE - HYBRID EKF DATA PROVENANCE'
        ['Exported: ' datestr(now,31)]
        ['Status: ' status]
        'All numbers come from the loaded Hybrid_EKF_Integrity_results.mat, never generated by this exporter.'
        sprintf('Epochs: evaluated=%d, valid=%d; burn-in=%d',m.nAssessed,m.nValid,settings.burnInEpochs)
        sprintf('GNSS sigma=%.9g m, TDOA difference sigma=%.9g m',settings.sigmaGNSS_m,settings.sigmaTDOA_m)
        sprintf('Nominal epoch risk allocation total=%.9g (horizontal+vertical)',settings.pHMI_per_epoch)
        sprintf('Per-update NIS false-alarm probability=%.9g',settings.pFalseAlarm)
        sprintf('Illustrative HAL=%.9g m, VAL=%.9g m',settings.HAL_m,settings.VAL_m)
        sprintf('TDOA reference metadata present=%d (not independently verified by presence alone)',referenceProvided)
        sprintf('Initial position from simulation truth=%d',truthInit)
        sprintf('NIS measurement block gating enabled=%d',gate)
        sprintf('Injected fault enabled=%d',faultEnabled)
        sprintf('Positive deterministic bias allowances enabled=%d',boundedMode)
        'Gaussian bounds are model-conditional. Gating introduces selection; no formal fault-risk allocation is implemented.'
        'A NIS event-level threshold is not an end-to-end false-alarm or missed-detection probability.'
        'Consecutive EKF epochs are correlated; sample counts are not independent rare-event trials.'
        'Figures: vector PDF + 600 dpi PNG, 17.8 cm width; no titles or captions embedded.'
        'Tables: booktabs environments, no vertical rules, captions in LaTeX.'
        'Outputs for Galileo LS / RAIM and standalone 5G WLS are NOT generated here.'
        'Reporting must state any unverified 5G reference mapping and truth-assisted initialization.'
        'To regenerate: S=load(''Hybrid_EKF_Integrity_results.mat''); Hybrid_EKF_Publication_Export(S.results,S.settings,''springer_export'');'
    };
    if isfield(settings,'Qmode')
        L{end+1}=['Process noise mode: ' char(settings.Qmode)];
    end
    if isfield(settings,'tdoaUnits')
        L{end+1}=['TDOA units: ' char(settings.tdoaUnits)];
    end
    if isfield(settings,'tdoaSign')
        L{end+1}=sprintf('TDOA sign convention: %.0f',settings.tdoaSign);
    end
    if isfield(settings,'correlatedTDOA')
        L{end+1}=sprintf('TDOA common-reference covariance model enabled: %d',logical(settings.correlatedTDOA));
    end
    if isfield(results,'sampleRates')
        sr=results.sampleRates;
        L{end+1}=sprintf('Sampling: receiver=%.9g Hz, Galileo=%.9g Hz, 5G period=%.9g s',...
            sr.receiver_Hz,sr.galileo_Hz,sr.nrPeriod_s);
    end
    if isfield(results,'sourceFiles')
        srcNames=fieldnames(results.sourceFiles);
        for i=1:numel(srcNames)
            v=results.sourceFiles.(srcNames{i});
            L{end+1}=['Source ' srcNames{i} ': ' char(v)]; %#ok<AGROW>
        end
    end
    writeTextLines(file,L);
end
function f=makeFigure(widthCm,heightCm)
    f=figure('Visible','off','Color','white','Units','centimeters', ...
             'Position',[2 2 widthCm heightCm], ...
             'PaperUnits','centimeters','PaperPosition',[0 0 widthCm heightCm]);
end
function formatAxes(f,showLegend)
    if nargin<2,showLegend=true;end
    ax=findall(f,'Type','axes');
    for k=1:numel(ax)
        set(ax(k),'FontName','Arial','FontSize',8,'LineWidth',0.65, ...
            'TickDir','out','Box','off');
        ax(k).XLabel.FontSize=9;
        ax(k).YLabel.FontSize=9;
        grid(ax(k),'on');
    end
    if showLegend
        for k=1:numel(ax)
            if numel(ax)==1
                legend(ax(k),'Location','best','Box','off','FontSize',8);
            end
        end
    end
end
function exportFigure(f,folder,basename)
    drawnow;
    exportgraphics(f,fullfile(folder,[basename '.pdf']), ...
                   'ContentType','vector','BackgroundColor','white');
    exportgraphics(f,fullfile(folder,[basename '.png']), ...
                   'Resolution',600,'BackgroundColor','white');
    % Optional figure source (.fig) can be regenerated if needed.
    close(f);
end
