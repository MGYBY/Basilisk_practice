function out = thixo_scan_FrVcrit_vs_T(a,Gamma,kappa,S0,Tvec,opts)
%THIXO_SCAN_FRVCRIT_VS_T  Long-wave neutral threshold FrVcrit as a function of T.
%
%   out = thixo_scan_FrVcrit_vs_T(a,Gamma,kappa,S0,Tvec,opts)
%
% Fixed physical/model parameters:
%   a      : thixotropic viscosity parameter, 0<a<=1
%   Gamma  : dimensionless shear-destruction parameter
%   kappa  : dimensionless structure diffusivity
%   S0     : slope parameter, S0=tan(theta)
%
% Scanned parameter:
%   Tvec   : vector of dimensionless relaxation times T
%
% Optional numerical/output settings in opts:
%   opts.N              Chebyshev order for long-wave linear BVPs [180]
%   opts.baseMeshN      initial bvp4c mesh size for base state [400]
%   opts.baseRelTol     base BVP relative tolerance [1e-5]
%   opts.baseAbsTol     base BVP absolute tolerance [1e-7]
%   opts.baseNMax       maximum base BVP mesh points [50000]
%   opts.continuation   kappa-continuation for base BVP [true]
%   opts.kappaStart     first continuation kappa [max(1e-2,kappa)]
%   opts.kappaSteps     number of continuation steps [] = automatic
%   opts.bvpSolver      'bvp4c' or 'bvp5c' ['bvp4c']
%   opts.diagnosticFrV  FrV used only to evaluate sigma2 diagnostic [4.8246]
%   opts.thresholdTol   tolerance for threshold classification [1e-10]
%   opts.outputPrefix   prefix for CSV and PNG outputs ['thixo_FrVcrit_vs_T']
%   opts.makePlot       true/false [true]
%   opts.writeCsv       true/false [true]
%   opts.verbose        true/false [true]
%
% Mathematical convention:
%   The complete kernel thixo_longwave_threshold_fixed.m computes
%
%       sigma(alpha) = -1i*c*alpha + Sigma2(FrV)*alpha^2 + O(alpha^3),
%       Sigma2(FrV) = Aq*FrV^2 + Bq.
%
%   Long waves are unstable when Sigma2>0.  A finite positive neutral
%   threshold exists only if
%
%       qcrit = FrVcrit^2 = -Bq/Aq > 0.
%
%   If qcrit<=0, the returned FrVcrit value is NaN and the status string
%   explains whether the branch is unstable for all FrV, stable for all FrV,
%   or numerically degenerate.
%
% Notes:
%   For fixed a,Gamma,kappa,S0, the base state is independent of T and FrV.
%   This scanner therefore solves the base-state BVP once, caches it, and
%   reuses it for all T values. This is faster and reduces numerical scatter
%   compared with calling the one-case code independently for every T.

    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    opts = local_defaults(opts,kappa);

    validate_inputs(a,Gamma,kappa,S0,Tvec);
    Tvec = Tvec(:);
    nT = numel(Tvec);

    if opts.verbose
        fprintf('\n--- Scanning long-wave critical FrV versus T ---\n');
        fprintf('Using kernel: %s\n', which('thixo_longwave_threshold_fixed'));
        fprintf('Fixed parameters: a=%g, Gamma=%g, kappa=%g, S0=%g\n', a,Gamma,kappa,S0);
        fprintf('T range: [%g, %g], nT=%d\n', min(Tvec), max(Tvec), nT);
    end

    FrVcrit = NaN(nT,1);
    ReVcrit = NaN(nT,1);
    qcrit = NaN(nT,1);
    Aq = NaN(nT,1);
    Bq = NaN(nT,1);
    phaseSpeed = NaN(nT,1);
    sigma2Diagnostic = NaN(nT,1);
    hasFiniteThreshold = false(nT,1);
    status = cell(nT,1);
    message = cell(nT,1);

    baseCache = [];

    for j = 1:nT
        T = Tvec(j);
        p = struct();
        p.a = a;
        p.T = T;
        p.Gamma = Gamma;
        p.kappa = kappa;
        p.S0 = S0;
        p.FrV = opts.diagnosticFrV;

        p.N = opts.N;
        p.baseMeshN = opts.baseMeshN;
        p.baseRelTol = opts.baseRelTol;
        p.baseAbsTol = opts.baseAbsTol;
        p.baseNMax = opts.baseNMax;
        p.continuation = opts.continuation;
        p.kappaStart = opts.kappaStart;
        p.kappaSteps = opts.kappaSteps;
        p.bvpSolver = opts.bvpSolver;
        p.mobilityFloor = opts.mobilityFloor;
        p.growthTol = opts.growthTol;
        p.newtonian = false;
        p.verbose = false;

        if ~isempty(baseCache)
            p.base = baseCache;
        end

        try
            kernel = thixo_longwave_threshold_fixed(p);
            if isempty(baseCache) && isfield(kernel,'base') && ~isempty(kernel.base)
                baseCache = kernel.base;
            end

            Aq(j) = kernel.A;
            Bq(j) = kernel.B;
            phaseSpeed(j) = kernel.c;
            sigma2Diagnostic(j) = kernel.sigma2;

            cls = classify_threshold(Aq(j),Bq(j),S0,opts.thresholdTol);
            hasFiniteThreshold(j) = cls.hasFiniteThreshold;
            FrVcrit(j) = cls.FrVcritNumeric;
            ReVcrit(j) = cls.ReVcritNumeric;
            qcrit(j) = cls.qcrit;
            status{j} = cls.status;
            message{j} = cls.message;

            if opts.verbose
                if cls.hasFiniteThreshold
                    fprintf('  %3d/%3d: T=%10.5g, FrVcrit=%12.6g, Aq=% .4e, Bq=% .4e, %s\n', ...
                        j,nT,T,FrVcrit(j),Aq(j),Bq(j),cls.status);
                else
                    fprintf('  %3d/%3d: T=%10.5g, FrVcrit=NaN, Aq=% .4e, Bq=% .4e, %s\n', ...
                        j,nT,T,Aq(j),Bq(j),cls.status);
                end
            end
        catch ME
            status{j} = 'error';
            message{j} = ME.message;
            if opts.verbose
                fprintf('  %3d/%3d: T=%10.5g, ERROR: %s\n', j,nT,T,ME.message);
            end
        end
    end

    Tcol = Tvec;
    ResultTable = table(Tcol,FrVcrit,ReVcrit,qcrit,Aq,Bq,phaseSpeed, ...
        sigma2Diagnostic,hasFiniteThreshold,status,message, ...
        'VariableNames',{'T','FrVcrit','ReVcrit','qcrit','Aq','Bq', ...
        'phaseSpeed_c','Sigma2_at_diagnosticFrV','hasFiniteThreshold','status','message'});

    out = struct();
    out.codeVersion = 'FRVCRIT_VS_T_SCAN_2026_05_10';
    out.input = struct('a',a,'Gamma',Gamma,'kappa',kappa,'S0',S0,'Tvec',Tvec);
    out.opts = opts;
    out.table = ResultTable;
    out.T = Tvec;
    out.FrVcrit = FrVcrit;
    out.ReVcrit = ReVcrit;
    out.qcrit = qcrit;
    out.Aq = Aq;
    out.Bq = Bq;
    out.phaseSpeed_c = phaseSpeed;
    out.sigma2Diagnostic = sigma2Diagnostic;
    out.hasFiniteThreshold = hasFiniteThreshold;
    out.status = status;
    out.message = message;
    out.baseCache = baseCache;

    if opts.writeCsv
        csvFile = [opts.outputPrefix '_scan.csv'];
        writetable(ResultTable,csvFile);
        out.csvFile = csvFile;
        if opts.verbose
            fprintf('\nWrote CSV: %s\n', csvFile);
        end
    else
        out.csvFile = '';
    end

    if opts.makePlot
        out.plotFile = make_threshold_plot(Tvec,FrVcrit,status,opts);
        out.coeffPlotFile = make_coeff_plot(Tvec,Aq,Bq,opts);
        if opts.verbose
            fprintf('Wrote plot: %s\n', out.plotFile);
            fprintf('Wrote coefficient plot: %s\n', out.coeffPlotFile);
        end
    else
        out.plotFile = '';
        out.coeffPlotFile = '';
    end

    if opts.verbose
        nFinite = sum(hasFiniteThreshold);
        nError = sum(strcmp(status,'error'));
        fprintf('\nScan completed: finite thresholds at %d/%d T values; errors at %d/%d T values.\n', ...
            nFinite,nT,nError,nT);
    end
end

% ----------------------------------------------------------------------- %
function opts = local_defaults(opts,kappa)
    d = struct();
    d.N = 180;
    d.baseMeshN = 400;
    d.baseRelTol = 1e-5;
    d.baseAbsTol = 1e-7;
    d.baseNMax = 50000;
    d.continuation = true;
    d.kappaStart = max(1e-2,kappa);
    d.kappaSteps = [];
    d.bvpSolver = 'bvp4c';
    d.mobilityFloor = 1e-14;
    d.growthTol = 1e-10;
    d.thresholdTol = 1e-10;
    d.diagnosticFrV = 4.8246;
    d.outputPrefix = 'thixo_FrVcrit_vs_T';
    d.makePlot = true;
    d.writeCsv = true;
    d.verbose = true;
    d.xScale = 'linear';   % 'linear' or 'log'
    names = fieldnames(d);
    for k = 1:numel(names)
        nm = names{k};
        if ~isfield(opts,nm) || isempty(opts.(nm))
            opts.(nm) = d.(nm);
        end
    end
end

% ----------------------------------------------------------------------- %
function validate_inputs(a,Gamma,kappa,S0,Tvec)
    if ~(isscalar(a) && isfinite(a) && a > 0 && a <= 1)
        error('Require scalar 0<a<=1.');
    end
    if ~(isscalar(Gamma) && isfinite(Gamma) && Gamma >= 0)
        error('Require scalar Gamma>=0.');
    end
    if ~(isscalar(kappa) && isfinite(kappa) && kappa > 0)
        error('Require scalar kappa>0.');
    end
    if ~(isscalar(S0) && isfinite(S0) && S0 > 0)
        error('Require scalar S0>0.');
    end
    if isempty(Tvec) || any(~isfinite(Tvec(:))) || any(Tvec(:) <= 0)
        error('Tvec must contain positive finite values.');
    end
end

% ----------------------------------------------------------------------- %
function cls = classify_threshold(Aq,Bq,S0,tol)
    scale = max([1,abs(Aq),abs(Bq)]);
    small = tol*scale;

    cls = struct();
    cls.hasFiniteThreshold = false;
    cls.FrVcritNumeric = NaN;
    cls.ReVcritNumeric = NaN;
    cls.qcrit = NaN;
    cls.status = 'undetermined';
    cls.message = '';

    if ~isfinite(Aq) || ~isfinite(Bq)
        cls.status = 'invalid_coefficients';
        cls.message = 'Aq or Bq is not finite.';
        return;
    end

    if abs(Aq) <= small
        if Bq > small
            cls.status = 'unstable_all_FrV_no_finite_threshold';
            cls.message = 'Sigma2 is positive and nearly independent of FrV.';
        elseif Bq < -small
            cls.status = 'stable_all_FrV_no_threshold';
            cls.message = 'Sigma2 is negative and nearly independent of FrV.';
        else
            cls.status = 'near_zero_all_FrV';
            cls.message = 'Sigma2 is numerically close to zero for all FrV.';
        end
        return;
    end

    qcrit = -Bq/Aq;
    cls.qcrit = qcrit;

    if isfinite(qcrit) && qcrit > 0
        cls.hasFiniteThreshold = true;
        cls.FrVcritNumeric = sqrt(qcrit);
        cls.ReVcritNumeric = qcrit/S0;
        if Aq > 0
            cls.status = 'finite_threshold_stable_below_unstable_above';
            cls.message = 'Finite neutral threshold: Sigma2<0 below and Sigma2>0 above.';
        else
            cls.status = 'finite_threshold_unstable_below_stable_above';
            cls.message = 'Finite neutral threshold: Sigma2>0 below and Sigma2<0 above.';
        end
    else
        if Aq > 0 && Bq >= -small
            cls.status = 'unstable_all_FrV_no_finite_threshold';
            cls.message = 'Sigma2(0)>=0 and Aq>0; no finite positive neutral FrV.';
        elseif Aq < 0 && Bq <= small
            cls.status = 'stable_all_FrV_no_threshold';
            cls.message = 'Sigma2(0)<=0 and Aq<0; no finite positive neutral FrV.';
        else
            cls.status = 'no_positive_root_roundoff_sensitive';
            cls.message = 'No positive qcrit found; classification may be roundoff-sensitive.';
        end
    end
end

% ----------------------------------------------------------------------- %
function plotFile = make_threshold_plot(Tvec,FrVcrit,status,opts)
    plotFile = [opts.outputPrefix '_FrVcrit_vs_T.png'];
    finiteMask = isfinite(FrVcrit);

    fig = figure('Visible','off');
    hold on; box on; grid on;

    if any(finiteMask)
        if strcmpi(opts.xScale,'log')
            semilogx(Tvec(finiteMask),FrVcrit(finiteMask),'o-','LineWidth',1.5,'MarkerSize',5);
        else
            plot(Tvec(finiteMask),FrVcrit(finiteMask),'o-','LineWidth',1.5,'MarkerSize',5);
        end
    else
        plot(Tvec,NaN(size(Tvec)),'o-','LineWidth',1.5,'MarkerSize',5);
        xlim([min(Tvec), max(Tvec)]);
        ylim([0,1]);
        text(mean([min(Tvec),max(Tvec)]),0.55, ...
            'No finite positive FrV_{crit} in this T range', ...
            'HorizontalAlignment','center','FontSize',12);
    end

    xlabel('relaxation time T');
    ylabel('critical Froude number FrV_{crit}');
    title('Long-wave neutral threshold FrV_{crit} versus T');

    noFiniteMask = ~finiteMask & ~strcmp(status,'error');
    if any(noFiniteMask) && any(finiteMask)
        yl = ylim;
        ymark = yl(1) + 0.03*(yl(2)-yl(1));
        plot(Tvec(noFiniteMask), ymark*ones(sum(noFiniteMask),1), 'x', 'MarkerSize',5);
        legend('finite threshold','no finite threshold','Location','best');
    elseif any(finiteMask)
        legend('finite threshold','Location','best');
    end

    saveas(fig,plotFile);
    close(fig);
end

% ----------------------------------------------------------------------- %
function coeffPlotFile = make_coeff_plot(Tvec,Aq,Bq,opts)
    coeffPlotFile = [opts.outputPrefix '_coefficients_vs_T.png'];
    fig = figure('Visible','off');

    subplot(2,1,1);
    plot(Tvec,Aq,'o-','LineWidth',1.3,'MarkerSize',4);
    grid on; box on;
    xlabel('relaxation time T');
    ylabel('A_q');
    title('\Sigma_2(FrV)=A_q FrV^2 + B_q');

    subplot(2,1,2);
    plot(Tvec,Bq,'o-','LineWidth',1.3,'MarkerSize',4);
    grid on; box on;
    xlabel('relaxation time T');
    ylabel('B_q');

    saveas(fig,coeffPlotFile);
    close(fig);
end
