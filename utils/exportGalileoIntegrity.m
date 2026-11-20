function summary = exportGalileoIntegrity(D, cfg)
%EXPORTGALILEOINTEGRITY Prepare publication outputs for Galileo RAIM.
%
%   summary = exportGalileoIntegrity(D, cfg)
%
% Mandatory D fields (all aligned to the same Galileo epoch):
%   t_s        [N x 1] time since start [s]
%   trueENU    [N x 3] ground-truth ENU position [m]
%   estENU     [N x 3] RAIM accepted/final position [m]
%   raimStatus [N x 1] status: PASS/EXCLUDED/FAULT_DETECTED/...
%   HPL        [N x 1] PRECOMPUTED horizontal protection level [m]
%
% Optional D fields:
%   nVisible   [N x 1] number of initially visible satellites
%   Tinitial   [N x 1] ORIGINAL global chi-square statistic
%   excludedSat[N x 1] index of excluded satellite (NaN if none)
%
% Mandatory assumptions for evaluation:
% - D.trueENU and D.estENU must use the SAME ENU reference origin.
% - D.trueENU must be aligned using actual timestamps or a verified
%   common epoch start and sample rate, never just matching array lengths.
% - Only the PASS branch is assigned a presently interpretable PL.
% - The resulting HMI frequencies are DESCRIPTIVE, NOT integrity guarantees.
%
% cfg fields (all optional; defaults below):
%   outDir, profile ('springer' or 'nature'), HAL_m [5 .5],
%   scenarioNames, sigma_m, Pfa, Pmd, alphaH, alphaV
%
% MATLAB: requires exportgraphics (R2020a or later), plus the Statistics
% and Machine Learning Toolbox for chi2inv if nVisible is supplied.
%
% Outputs: Tables/*.csv and *.tex; Figures/*.pdf, *.eps, *.png, *.fig.

    if nargin < 2 || isempty(cfg), cfg = struct(); end
    cfg = defaults(cfg);

    needed = {'t_s','trueENU','estENU','raimStatus','HPL'};
    for i = 1:numel(needed)
        assert(isfield(D,needed{i}), 'Missing D.%s', needed{i});
    end

    t = double(D.t_s(:));
    N = numel(t);
    assert(N > 0 && all(isfinite(t)) && all(diff(t) > 0), ...
        'D.t_s must be a nonempty, strictly increasing time vector.');
    assert(isequal(size(D.trueENU),[N 3]),'D.trueENU must be N-by-3.');
    assert(isequal(size(D.estENU),[N 3]),'D.estENU must be N-by-3.');

    stat = string(D.raimStatus(:));
    hpl = double(D.HPL(:));
    assert(numel(stat)==N && numel(hpl)==N, 'Inconsistent epoch lengths.');
    assert(numel(cfg.HAL_m)==2 && numel(cfg.scenarioNames)==2, ...
        'This exporter expects exactly two application scenarios.');
    assert(all(cfg.HAL_m>0), 'HAL must be positive.');

    satN = nan(N,1);
    if isfield(D,'nVisible')
        satN = double(D.nVisible(:));
        assert(numel(satN)==N, 'D.nVisible must have N elements.');
    end
    Tinitial = nan(N,1);
    if isfield(D,'Tinitial')
        Tinitial = double(D.Tinitial(:));
        assert(numel(Tinitial)==N, 'D.Tinitial must have N elements.');
    end
    excludedSat = nan(N,1);
    if isfield(D,'excludedSat')
        excludedSat = double(D.excludedSat(:));
        assert(numel(excludedSat)==N, 'D.excludedSat must have N elements.');
    end

    % Explicitly use pre-FDE test degrees of freedom for the original test.
    gammaInitial = nan(N,1);
    if any(isfinite(satN))
        testable = isfinite(satN) & satN >= 5;
        gammaInitial(testable) = chi2inv(1-cfg.Pfa,satN(testable)-4);
    end

    truthValid = all(isfinite(D.trueENU(:,1:2)),2);
    estimateValid = all(isfinite(D.estENU(:,1:2)),2);
    hpe = nan(N,1);
    validError = truthValid & estimateValid;
    errH = D.estENU(validError,1:2)-D.trueENU(validError,1:2);
    hpe(validError) = sqrt(sum(errH.^2,2));

    monitorPASS = stat == "PASS" & estimateValid & isfinite(hpl) & hpl>=0;
    nSc = numel(cfg.HAL_m);
    usable = false(N,nSc);
    HMI = false(N,nSc);
    plViolation = false(N,nSc);

    nAvail = zeros(nSc,1);
    nHMI = zeros(nSc,1);
    nPLviol = zeros(nSc,1);
    nAcceptedNoTruth = zeros(nSc,1);
    medHPL = nan(nSc,1);
    p95HPL = nan(nSc,1);
    medHPE = nan(nSc,1);
    p95HPE = nan(nSc,1);

    for s = 1:nSc
        usable(:,s) = monitorPASS & hpl <= cfg.HAL_m(s);
        HMI(:,s) = usable(:,s) & truthValid & hpe > cfg.HAL_m(s);
        plViolation(:,s) = usable(:,s) & truthValid & hpe > hpl;
        nAvail(s) = nnz(usable(:,s));
        nHMI(s) = nnz(HMI(:,s));
        nPLviol(s) = nnz(plViolation(:,s));
        nAcceptedNoTruth(s) = nnz(usable(:,s) & ~truthValid);
        medHPL(s) = medianOrNaN(hpl(usable(:,s)),50);
        p95HPL(s) = medianOrNaN(hpl(usable(:,s)),95);
        medHPE(s) = medianOrNaN(hpe(usable(:,s) & truthValid),50);
        p95HPE(s) = medianOrNaN(hpe(usable(:,s) & truthValid),95);
    end

    nTruth = nnz(truthValid);
    observedHMIperTruthEpoch = nan(nSc,1);
    if nTruth > 0
        observedHMIperTruthEpoch = nHMI / nTruth;
    end

    summary = table(cfg.scenarioNames(:),cfg.HAL_m(:), ...
        repmat(N,nSc,1), repmat(nTruth,nSc,1), nAvail, ...
        100*nAvail/N, nHMI, observedHMIperTruthEpoch, ...
        nPLviol, nAcceptedNoTruth, medHPL, p95HPL, ...
        medHPE, p95HPE, ...
        'VariableNames', {'Scenario','HAL_m','N_epochs','N_truth', ...
        'N_available','Availability_pct','N_HMI_observed', ...
        'Observed_HMI_per_truth_epoch','N_PL_exceed', ...
        'N_available_missing_truth','Median_HPL_m', 'P95_HPL_m', ...
        'Median_HPE_m','P95_HPE_m'});

    tableDir = fullfile(cfg.outDir,'Tables');
    figDir = fullfile(cfg.outDir,'Figures');
    if ~exist(tableDir,'dir'), mkdir(tableDir); end
    if ~exist(figDir,'dir'), mkdir(figDir); end

    % Raw, traceable scientific outputs.
    writetable(summary,fullfile(tableDir,'Table1_IntegritySummary.csv'));
    epochTable = table((1:N).',t,stat,satN,excludedSat,Tinitial, ...
        gammaInitial,hpe,hpl,usable(:,1),usable(:,2), ...
        HMI(:,1),HMI(:,2),plViolation(:,1),plViolation(:,2), ...
        'VariableNames',{'Epoch','Time_s','RAIM_status','N_visible', ...
        'Excluded_sat','T_initial','Gamma_initial','HPE_m','HPL_m', ...
        'Available_A','Available_B','ObservedHMI_A','ObservedHMI_B', ...
        'PLViolation_A','PLViolation_B'});
    writetable(epochTable,fullfile(tableDir,'GalileoEpochResults.csv'));

    % RAIM status distribution.
    statuses = ["PASS";"EXCLUDED";"FAULT_DETECTED"; ...
                "AMBIGUOUS";"UNAVAILABLE"];
    statusCount = arrayfun(@(k)nnz(stat==statuses(k)),(1:numel(statuses)).');
    statusTable = table(statuses,statusCount,100*statusCount/N,...
        'VariableNames',{'Status','N_epochs','Percentage'});
    writetable(statusTable,fullfile(tableDir,'Table2_RaimStatus.csv'));

    writeSummaryTex(fullfile(tableDir,'Table1_IntegritySummary.tex'),summary,nTruth);
    writeStatusTex(fullfile(tableDir,'Table2_RaimStatus.tex'),statusTable);
    writeConfigTex(fullfile(tableDir,'Table3_Configuration.tex'),cfg,N);
    writeCaptions(fullfile(cfg.outDir,'FigureCaptions.txt'));

    sty = graphicsStyle(cfg.profile);

    % Figure 1: pre-exclusion global RAIM T and threshold.
    validTest = isfinite(Tinitial) & isfinite(gammaInitial);
    if any(validTest)
        f = newFigure(sty,sty.doubleCm,7.5);
        ax = axes('Parent',f); hold(ax,'on'); grid(ax,'on');
        plot(ax,t(validTest),Tinitial(validTest),'-', ...
            'Color',sty.blue,'LineWidth',0.9,'DisplayName','Global RAIM statistic');
        plot(ax,t(validTest),gammaInitial(validTest),'--', ...
            'Color',sty.black,'LineWidth',1.0,'DisplayName','Detection threshold');
        triggered = validTest & Tinitial > gammaInitial;
        if any(triggered)
            scatter(ax,t(triggered),Tinitial(triggered),22,sty.orange,'x', ...
                'LineWidth',1.1,'DisplayName','Global test rejection');
        end
        xlabel(ax,'Time (s)'); ylabel(ax,'Normalized residual statistic');
        legend(ax,'Location','best');
        saveFigure(f,figDir,'Fig1_RaimTest',sty);
    else
        warning('Skipping Fig1: provide D.Tinitial and D.nVisible.');
    end

    % Figure 2: Stanford-style HPE/HPL, one panel per ALERT LIMIT.
    f = newFigure(sty,sty.doubleCm,8.8);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    maxVal = max([cfg.HAL_m(:);hpe(isfinite(hpe));hpl(isfinite(hpl));1]);
    xLow = 1e-3; xHigh = max(2*maxVal,10);
    for s = 1:nSc
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
        m = monitorPASS & truthValid;
        unavail = m & ~usable(:,s);
        safe = usable(:,s) & truthValid & hpe <= cfg.HAL_m(s);
        harm = HMI(:,s);
        x = max(hpe,xLow); y = max(hpl,xLow);
        if any(unavail)
            scatter(ax,x(unavail),y(unavail),16,sty.gray,'o', ...
                'DisplayName','Not available');
        end
        if any(safe)
            scatter(ax,x(safe),y(safe),18,sty.blue,'o','filled', ...
                'DisplayName','Available, HPE <= AL');
        end
        if any(harm)
            scatter(ax,x(harm),y(harm),32,sty.orange,'x', ...
                'LineWidth',1.2,'DisplayName','Observed HMI');
        end
        plot(ax,[xLow xHigh],[xLow xHigh],':','Color',sty.black, ...
            'LineWidth',0.8,'HandleVisibility','off');
        xline(ax,cfg.HAL_m(s),'--','Color',sty.black, ...
            'HandleVisibility','off');
        yline(ax,cfg.HAL_m(s),'--','Color',sty.black, ...
            'HandleVisibility','off');
        set(ax,'XScale','log','YScale','log', ...
            'XLim',[xLow xHigh],'YLim',[xLow xHigh]);
        axis(ax,'square');
        xlabel(ax,'HPE (m)'); ylabel(ax,'HPL (m)');
        title(ax,sprintf('(%c) AL = %.1f m','a'+s-1,cfg.HAL_m(s)), ...
            'FontWeight','normal');
        if s == 1 && (~isempty(findobj(ax,'Type','Scatter')))
            legend(ax,'Location','northwest');
        end
    end
    saveFigure(f,figDir,'Fig2_Stanford',sty);

    % Figure 3: availability of PASS solutions over hypothetical HAL.
    f = newFigure(sty,sty.singleCm,7.0);
    ax = axes('Parent',f); hold(ax,'on'); grid(ax,'on');
    maxH = hpl(monitorPASS);
    if isempty(maxH), upperHAL = 100; else, upperHAL=max(100,max(maxH)*1.2); end
    HALgrid = logspace(-1,log10(upperHAL),240);
    availCurve = arrayfun(@(a)100*mean(monitorPASS & hpl<=a),HALgrid);
    semilogx(ax,HALgrid,availCurve,'-','Color',sty.blue, ...
        'LineWidth',1.25,'DisplayName','Galileo RAIM PASS');
    for s = 1:nSc
        plot(ax,cfg.HAL_m(s),100*nAvail(s)/N,'o', ...
            'MarkerSize',5,'MarkerFaceColor',sty.markers(s,:), ...
            'MarkerEdgeColor',sty.markers(s,:), ...
            'DisplayName',sprintf('Scenario %c','A'+s-1));
        xline(ax,cfg.HAL_m(s),':','Color',sty.markers(s,:), ...
            'HandleVisibility','off');
    end
    xlabel(ax,'Horizontal alert limit (m)');
    ylabel(ax,'Availability (%)'); ylim(ax,[0 100]);
    legend(ax,'Location','southeast');
    saveFigure(f,figDir,'Fig3_Availability',sty);

    % Figure 4: time trace of observed position error and computed HPL.
    f = newFigure(sty,sty.doubleCm,7.5);
    ax = axes('Parent',f); hold(ax,'on'); grid(ax,'on');
    xx = hpe; yy = hpl;
    xx(~monitorPASS) = NaN; yy(~monitorPASS) = NaN;
    plot(ax,t,xx,'-','LineWidth',1.0,'Color',sty.orange, ...
        'DisplayName','Observed HPE');
    plot(ax,t,yy,'-','LineWidth',1.0,'Color',sty.blue, ...
        'DisplayName','Computed HPL');
    for s = 1:nSc
        yline(ax,cfg.HAL_m(s),'--','Color',sty.markers(s,:), ...
            'LineWidth',0.8,'DisplayName',sprintf('AL, scenario %c','A'+s-1));
    end
    xlabel(ax,'Time (s)'); ylabel(ax,'Horizontal error / bound (m)');
    legend(ax,'Location','best');
    saveFigure(f,figDir,'Fig4_PositionError',sty);

    % MATLAB binary reproducibility snapshot: data and presentation settings.
    save(fullfile(cfg.outDir,'IntegrityExportData.mat'), ...
        'summary','epochTable','statusTable','cfg','-v7');

    fprintf('Publication export written to: %s\n',cfg.outDir);
    disp(summary);
    if any(nAcceptedNoTruth>0)
        warning('Some available epochs lack truth. Observed HMI counts are incomplete.');
    end
    if nTruth < N
        warning('Only %d/%d epochs have ground truth. Report evaluation coverage.',nTruth,N);
    end
    fprintf(['WARNING: Observed HMI counts do not demonstrate any ' ...
        'integrity risk target; validate PL models and use fault Monte Carlo.\n']);
end

function cfg = defaults(cfg)
    if ~isfield(cfg,'outDir'), cfg.outDir='GalileoPublication'; end
    if ~isfield(cfg,'profile'), cfg.profile='springer'; end
    if ~isfield(cfg,'HAL_m'), cfg.HAL_m=[5;0.5]; end
    if ~isfield(cfg,'scenarioNames')
        cfg.scenarioNames=["A - Urban ITS";"B - Strict positioning"];
    end
    if ~isfield(cfg,'sigma_m'), cfg.sigma_m=4.6; end
    if ~isfield(cfg,'Pfa'), cfg.Pfa=1e-3; end
    if ~isfield(cfg,'Pmd'), cfg.Pmd=1e-3; end
    if ~isfield(cfg,'alphaH'), cfg.alphaH=1e-5; end
    if ~isfield(cfg,'alphaV'), cfg.alphaV=1e-5; end
    cfg.HAL_m=cfg.HAL_m(:);
    cfg.scenarioNames=string(cfg.scenarioNames(:));
    assert(cfg.Pfa>0 && cfg.Pfa<1);
    assert(cfg.Pmd>0 && cfg.Pmd<1);
    assert(any(strcmpi(cfg.profile,{'springer','nature'})), ...
        'cfg.profile must be springer or nature.');
end

function style = graphicsStyle(profile)
    style.blue   = [0 0.4470 0.7410];
    style.orange = [0.8500 0.3250 0.0980];
    style.black  = [0.20 0.20 0.20];
    style.gray   = [0.55 0.55 0.55];
    style.markers = [style.blue;style.orange];
    if strcmpi(profile,'nature')
        style.singleCm=8.9;
        style.doubleCm=18.3;
        style.fontPt=6.5;
    else
        style.singleCm=8.4;
        style.doubleCm=17.4;
        style.fontPt=8;
    end
end

function f = newFigure(style,widthCm,heightCm)
    f=figure('Visible','off','Color','w','Units','centimeters', ...
        'Position',[2 2 widthCm heightCm], ...
        'PaperPositionMode','auto','Renderer','painters');
    set(f,'DefaultAxesFontName','Arial', ...
        'DefaultTextFontName','Arial', ...
        'DefaultAxesFontSize',style.fontPt, ...
        'DefaultTextFontSize',style.fontPt);
end

function saveFigure(f,figDir,base,style)
    ax=findall(f,'Type','axes');
    set(ax,'FontName','Arial','FontSize',style.fontPt,'LineWidth',0.6, ...
        'Box','off','TickDir','out');
    exportgraphics(f,fullfile(figDir,[base '.pdf']), ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(f,fullfile(figDir,[base '.png']), ...
        'Resolution',600,'BackgroundColor','white');
    savefig(f,fullfile(figDir,[base '.fig']));
    try
        print(f,fullfile(figDir,base),'-depsc','-painters');
    catch ME
        warning('EPS export failed for %s: %s',base,ME.message);
    end
    close(f);
end

function x = medianOrNaN(a,p)
    a=a(isfinite(a));
    if isempty(a), x=NaN; else, x=prctile(a,p); end
end

function writeSummaryTex(file,tab,nTruth)
    fid = fopen(file,'w'); assert(fid>0,'Cannot create LaTeX table.');
    cleaner=onCleanup(@()fclose(fid));
    fprintf(fid,'%% Requires \\usepackage{booktabs} in manuscript preamble.\n');
    fprintf(fid,'\\begin{table*}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Galileo RAIM integrity evaluation under two horizontal alert limits.}\n');
    fprintf(fid,'\\label{tab:galileo-integrity}\n');
    fprintf(fid,'\\begin{tabular}{lrrrrrr}\n\\toprule\n');
    fprintf(fid,'Scenario & AL (m) & Availability (\\%%) & HMI (n) & PL exc. (n) & Median HPL (m) & P95 HPL (m) \\\\\n');
    fprintf(fid,'\\midrule\n');
    for i=1:height(tab)
        fprintf(fid,'%s & %.2f & %.2f & %d & %d & %s & %s \\\\\n', ...
            char(tab.Scenario(i)),tab.HAL_m(i),tab.Availability_pct(i), ...
            tab.N_HMI_observed(i),tab.N_PL_exceed(i), ...
            fmt(tab.Median_HPL_m(i)),fmt(tab.P95_HPL_m(i)));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
    fprintf(fid,'\\par\\smallskip\\footnotesize{');
    fprintf(fid,['Availability: accepted PASS epochs with HPL $\\leq$ AL divided by all epochs. ' ...
        'HMI: observed epochs with accepted position and HPE $>$ AL. ' ...
        'PL exc.: accepted epochs with HPE $>$ HPL. ' ...
        'Observed HMI counts are not certified integrity risks. ' ...
        'Ground-truth epochs: %d.}\n'],nTruth);
    fprintf(fid,'\\end{table*}\n');
end

function writeStatusTex(file,tab)
    fid=fopen(file,'w'); assert(fid>0); cleaner=onCleanup(@()fclose(fid));
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Distribution of Galileo RAIM-FDE monitoring states.}\n');
    fprintf(fid,'\\label{tab:raim-status}\n');
    fprintf(fid,'\\begin{tabular}{lrr}\n\\toprule\n');
    fprintf(fid,'Status & Epochs & Share (\\%%) \\\\\n\\midrule\n');
    for i=1:height(tab)
        fprintf(fid,'%s & %d & %.2f \\\\\n', ...
            strrep(char(tab.Status(i)),'_','\_'), ...
            tab.N_epochs(i),tab.Percentage(i));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

function writeConfigTex(file,cfg,N)
    fid=fopen(file,'w'); assert(fid>0); cleaner=onCleanup(@()fclose(fid));
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Illustrative Galileo RAIM monitoring configuration.}\n');
    fprintf(fid,'\\label{tab:galileo-config}\n');
    fprintf(fid,'\\begin{tabular}{lr}\n\\toprule\n');
    fprintf(fid,'Parameter & Value \\\\\n\\midrule\n');
    fprintf(fid,'Number of Galileo epochs & %d \\\\\n',N);
    fprintf(fid,'Common pseudorange $\\sigma$ (m) & %.3g \\\\\n',cfg.sigma_m);
    fprintf(fid,'Global test $P_{FA}$ & %.3g \\\\\n',cfg.Pfa);
    fprintf(fid,'Fault test $P_{MD}$ & %.3g \\\\\n',cfg.Pmd);
    fprintf(fid,'Nominal tail $\\alpha_H$ & %.3g \\\\\n',cfg.alphaH);
    fprintf(fid,'Nominal tail $\\alpha_V$ & %.3g \\\\\n',cfg.alphaV);
    fprintf(fid,'Scenario A AL (m) & %.3g \\\\\n',cfg.HAL_m(1));
    fprintf(fid,'Scenario B AL (m) & %.3g \\\\\n',cfg.HAL_m(2));
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

function writeCaptions(file)
    fid=fopen(file,'w'); assert(fid>0); cleaner=onCleanup(@()fclose(fid));
    fprintf(fid,['Fig. 1. Global Galileo RAIM chi-square statistic and the ' ...
        'epoch-dependent detection threshold before fault exclusion. ' ...
        'Crosses indicate rejected global tests.\n\n']);
    fprintf(fid,['Fig. 2. Stanford-style diagrams for (a) urban ITS (AL=5 m) ' ...
        'and (b) the strict horizontal AL sensitivity analysis (AL=0.5 m). ' ...
        'Shown points have valid ground truth, an estimated position, and a ' ...
        'computed PL from the global RAIM PASS branch. Dashed lines indicate ' ...
        'the AL; the dotted diagonal denotes HPE=HPL. Zero values are plotted ' ...
        'at 0.001 m solely for log-scale display.\n\n']);
    fprintf(fid,['Fig. 3. Empirical availability of Galileo global-RAIM PASS ' ...
        'solutions as a function of horizontal AL for an unchanged computed ' ...
        'PL. Markers show the two application scenarios.\n\n']);
    fprintf(fid,['Fig. 4. Time history of horizontal position error (HPE), ' ...
        'the computed protection level (HPL), and horizontal AL values. ' ...
        'Only epochs with a global RAIM PASS decision are shown.\n']);
end

function s=fmt(x)
    if ~isfinite(x), s='--'; else, s=sprintf('%.2f',x); end
end
