%RUN_THIXO_SCAN_FRVCRIT_VS_T
% Demonstration: long-wave critical FrV versus T for fixed thixotropic
% parameters.
%
% This script scans T from 0.1 to 2 for
%   a      = 0.2
%   Gamma  = 8
%   kappa  = 1e-4
%   S0     = 0.05
%
% It writes:
%   demo_FrVcrit_vs_T_scan.csv
%   demo_FrVcrit_vs_T_FrVcrit_vs_T.png
%   demo_FrVcrit_vs_T_coefficients_vs_T.png

clear; clc;

fprintf('Running FrVcrit-versus-T scan.\n');
fprintf('Scanner resolved by MATLAB:\n  %s\n', which('thixo_scan_FrVcrit_vs_T'));
fprintf('Kernel resolved by MATLAB:\n  %s\n\n', which('thixo_longwave_threshold_fixed'));

% Fixed case definition requested by the user.
a      = 0.2;
Gamma  = 8.0;
kappa  = 1.0e-4;
S0     = 0.05;

% Requested T range.
Tmin = 0.1;
Tmax = 2.0;
nT   = 50;
Tvec = linspace(Tmin,Tmax,nT);

% Numerical options. For kappa=1e-4 the base structure layer is thin, so
% N=180 is a reasonable first production setting. Increase to N=220 for a
% convergence audit after the first run.
opts = struct();
opts.N = 180;
opts.baseMeshN = 400;
opts.baseRelTol = 1e-5;
opts.baseAbsTol = 1e-7;
opts.baseNMax = 50000;
opts.continuation = true;
opts.kappaStart = 1e-2;
opts.kappaSteps = [];
opts.bvpSolver = 'bvp4c';
opts.outputPrefix = 'demo_FrVcrit_vs_T';
opts.makePlot = true;
opts.writeCsv = true;
opts.verbose = true;
opts.diagnosticFrV = 4.8246;

out = thixo_scan_FrVcrit_vs_T(a,Gamma,kappa,S0,Tvec,opts);

fprintf('\nSummary:\n');
fprintf('  finite thresholds: %d/%d T values\n', sum(out.hasFiniteThreshold), numel(out.T));
if any(out.hasFiniteThreshold)
    fprintf('  min finite FrVcrit = %.12g\n', min(out.FrVcrit(out.hasFiniteThreshold)));
    fprintf('  max finite FrVcrit = %.12g\n', max(out.FrVcrit(out.hasFiniteThreshold)));
else
    fprintf('  no finite positive FrVcrit values were found in this T range.\n');
end
fprintf('  CSV:  %s\n', out.csvFile);
fprintf('  plot: %s\n', out.plotFile);
fprintf('  coefficient plot: %s\n', out.coeffPlotFile);
