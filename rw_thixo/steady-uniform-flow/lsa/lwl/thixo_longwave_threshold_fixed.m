function out = thixo_longwave_threshold_fixed(p)
%THIXO_LONGWAVE_THRESHOLD_FIXED Corrected alpha->0 stability coefficient.
%
%   out = thixo_longwave_threshold_fixed(p)
%
% Computes the long-wave temporal dispersion expansion
%
%       sigma(alpha) = -1i*c*alpha + sigma2*alpha^2 + O(alpha^3),
%
% for the dimensionless thixotropic falling-film equations used in this
% project, without surface tension.  The calculation avoids the finite-alpha
% Orr--Sommerfeld spectrum and solves only one nonlinear base-state BVP plus
% linear wall-normal long-wave BVPs.
%
% Required physical input fields:
%       p.a       viscosity-contrast parameter, 0<a<=1
%       p.FrV     Froude number based on the project velocity scale V
%       p.T       dimensionless structural relaxation time
%       p.Gamma   dimensionless shear-destruction parameter
%       p.kappa   dimensionless structural diffusivity
%       p.S0      slope parameter, S0=tan(theta)
%
% Main numerical fields, all optional:
%       p.N              Chebyshev order for long-wave linear BVPs [180]
%       p.baseMeshN      initial bvp4c mesh size [400]
%       p.baseRelTol     base bvp relative tolerance [1e-5]
%       p.baseAbsTol     base bvp absolute tolerance [1e-7]
%       p.baseNMax       maximum base bvp mesh points [50000]
%       p.continuation   continue from larger kappa to target kappa [true]
%       p.kappaStart     first continuation kappa [max(1e-2,p.kappa)]
%       p.kappaSteps     number of continuation values; [] => automatic
%       p.bvpSolver      'bvp4c' or 'bvp5c' if available ['bvp4c']
%       p.verbose        print progress [true]
%       p.newtonian      analytic Newtonian validation mode [false]
%
% Important corrections relative to the earlier version:
%   1. Continuation now reuses the previous bvp4c solution through deval,
%      not through bvpinit(sol.x,sol.y), which is invalid when sol.y is a
%      matrix of values at mesh points.
%   2. The equations are written in terms of mobility A(lambda)=1/M(lambda),
%      avoiding the nearly singular viscosity M near lambda=1.
%   3. The O(alpha) correction uses the pressure-eliminated/vorticity form.
%      This recovers the exact Newtonian long-wave threshold
%          sigma2_N=(2*FrV^2-5)/(15*S0),  FrV_c=sqrt(5/2).
%
% Sign convention:
%       disturbance ~ exp(i*alpha*x + sigma*t)
%       long waves are unstable if sigma2 > 0.

    p = set_defaults(p);
    validate_params(p);

    if p.newtonian
        out = newtonian_result(p);
        return;
    end

    qFr = p.FrV^2;
    if p.verbose
        fprintf('\n--- Corrected thixotropic long-wave threshold calculation ---\n');
        fprintf('a=%g, FrV=%g, T=%g, Gamma=%g, kappa=%g, S0=%g\n', ...
            p.a, p.FrV, p.T, p.Gamma, p.kappa, p.S0);
        fprintf('Derived: beta=S0/FrV^2=%g, G=1/FrV^2=%g, rlambda=S0/T=%g, ReV=FrV^2/S0=%g\n', ...
            p.S0/qFr, 1/qFr, p.S0/p.T, qFr/p.S0);
    end

    if isfield(p,'base') && ~isempty(p.base)
        base = p.base;
        if ~isfield(base,'N') || base.N ~= p.N
            base = build_base_on_cheb_grid(p, base.sol);
        end
    else
        sol = solve_base_bvp(p);
        base = build_base_on_cheb_grid(p, sol);
    end

    lead = solve_leading_problem(p, base);

    % sigma2 is affine in q=FrV^2 for fixed (a,T,Gamma,kappa,S0):
    %     sigma2(q) = A*q + B.
    corr0 = solve_correction_problem(p, base, lead, 0.0);
    corr1 = solve_correction_problem(p, base, lead, 1.0);
    B = corr0.sigma2;
    Acoef = corr1.sigma2 - corr0.sigma2;
    sigma2 = Acoef*qFr + B;

    if Acoef ~= 0 && -B/Acoef > 0
        FrVcrit = sqrt(-B/Acoef);
        ReVcrit = FrVcrit^2/p.S0;
    else
        FrVcrit = NaN;
        ReVcrit = NaN;
    end

    if sigma2 > p.growthTol
        status = 'unstable';
    elseif sigma2 < -p.growthTol
        status = 'stable';
    else
        status = 'neutral/undetermined';
    end

    out = struct();
    out.params = p;
    out.beta = p.S0/qFr;
    out.G = 1/qFr;
    out.rlambda = p.S0/p.T;
    out.Tbar = p.T/p.S0;
    out.ReV = qFr/p.S0;
    out.c = lead.c;
    out.A = Acoef;
    out.B = B;
    out.sigma2 = sigma2;
    out.FrVcrit = FrVcrit;
    out.ReVcrit = ReVcrit;
    out.status = status;
    out.base = base;
    out.leading = lead;
    out.correction_q0 = corr0;
    out.correction_q1 = corr1;
    out.newtonian_sigma2_formula = (2*p.FrV^2 - 5)/(15*p.S0);

    if p.verbose
        fprintf('\nLong-wave result:\n');
        fprintf('  c                 = %.16g\n', out.c);
        fprintf('  sigma2(FrV)       = %.16g\n', out.sigma2);
        fprintf('  sigma2(q)=A*q+B   = (%.16g)*FrV^2 + (%.16g)\n', out.A, out.B);
        fprintf('  status            = %s\n', out.status);
        if ~isnan(out.FrVcrit)
            fprintf('  FrVcrit           = %.16g\n', out.FrVcrit);
            fprintf('  ReVcrit           = %.16g\n', out.ReVcrit);
        else
            fprintf('  FrVcrit           = NaN  (no positive affine neutral crossing)\n');
        end
    end
end

% ======================================================================= %
function p = set_defaults(p)
    defaults = struct();
    defaults.a = 0.2;
    defaults.FrV = 4.8246;
    defaults.T = 100.0;
    defaults.Gamma = 8.0;
    defaults.kappa = 1e-4;
    defaults.S0 = 0.05;

    defaults.N = 180;
    defaults.baseMeshN = 400;
    defaults.baseRelTol = 1e-5;
    defaults.baseAbsTol = 1e-7;
    defaults.baseNMax = 50000;
    defaults.continuation = true;
    defaults.kappaStart = [];
    defaults.kappaSteps = [];
    defaults.bvpSolver = 'bvp4c';
    defaults.mobilityFloor = 1e-14;   % only for diagnostic M=1/A output
    defaults.newtonian = false;
    defaults.verbose = true;
    defaults.growthTol = 1e-10;

    names = fieldnames(defaults);
    for k = 1:numel(names)
        nm = names{k};
        if ~isfield(p,nm) || isempty(p.(nm))
            p.(nm) = defaults.(nm);
        end
    end
    if isempty(p.kappaStart)
        p.kappaStart = max(1e-2, p.kappa);
    end
end

% ======================================================================= %
function validate_params(p)
    if p.a <= 0 || p.a > 1
        error('Require 0 < p.a <= 1.');
    end
    if p.FrV <= 0 || p.T <= 0 || p.Gamma < 0 || p.kappa <= 0 || p.S0 <= 0
        error('Require FrV>0, T>0, Gamma>=0, kappa>0 and S0>0.');
    end
    if p.N < 20
        error('p.N is too small. Use p.N >= 80 for serious calculations.');
    end
end

% ======================================================================= %
function out = newtonian_result(p)
    qFr = p.FrV^2;
    Acoef = 2/(15*p.S0);
    B = -1/(3*p.S0);
    sigma2 = Acoef*qFr + B;

    out = struct();
    out.params = p;
    out.beta = p.S0/qFr;
    out.G = 1/qFr;
    out.rlambda = p.S0/p.T;
    out.Tbar = p.T/p.S0;
    out.ReV = qFr/p.S0;
    out.c = 1.0;
    out.A = Acoef;
    out.B = B;
    out.sigma2 = sigma2;
    out.FrVcrit = sqrt(5/2);
    out.ReVcrit = out.FrVcrit^2/p.S0;
    out.status = ternary_status(sigma2,p.growthTol);
    out.base = [];
    out.leading = [];
    out.correction_q0 = [];
    out.correction_q1 = [];
    out.newtonian_sigma2_formula = sigma2;

    if p.verbose
        fprintf('\nNewtonian validation mode:\n');
        fprintf('  c = 1\n');
        fprintf('  sigma2 = (2*FrV^2 - 5)/(15*S0) = %.16g\n', sigma2);
        fprintf('  FrVcrit = sqrt(5/2) = %.16g\n', out.FrVcrit);
        fprintf('  ReVcrit = FrVcrit^2/S0 = %.16g\n', out.ReVcrit);
    end
end

function s = ternary_status(x,tol)
    if x > tol
        s = 'unstable';
    elseif x < -tol
        s = 'stable';
    else
        s = 'neutral/undetermined';
    end
end

% ======================================================================= %
function sol = solve_base_bvp(p)
    if p.verbose
        fprintf('\nSolving base state for Lambda(z) using mobility A=1/M...\n');
        fprintf('  target kappa = %.4e\n', p.kappa);
    end

    opts = bvpset('RelTol',p.baseRelTol, ...
                  'AbsTol',p.baseAbsTol, ...
                  'NMax',p.baseNMax, ...
                  'Stats','off');

    if p.continuation && p.kappa < p.kappaStart
        kappas = continuation_sequence(p.kappaStart, p.kappa, p.kappaSteps);
    else
        kappas = p.kappa;
    end

    solverName = lower(string(p.bvpSolver));
    useBvp5c = solverName == "bvp5c" && exist('bvp5c','file') == 2;

    sol = [];
    for j = 1:numel(kappas)
        pk = p;
        pk.kappa = kappas(j);
        zmesh = clustered_mesh(p.baseMeshN, pk.kappa);

        if isempty(sol)
            solinit = bvpinit(zmesh, @(z) base_initial_guess(z, pk));
        else
            % Correct continuation usage: use deval(sol,z), which returns a
            % vector at each requested z.  Passing sol.y directly to bvpinit
            % is invalid because sol.y is a matrix on the old adaptive mesh.
            solinit = bvpinit(zmesh, @(z) deval(sol,z));
        end

        if p.verbose
            fprintf('  base continuation step %2d/%2d: kappa = %.4e ... ', ...
                j, numel(kappas), pk.kappa);
        end

        if useBvp5c
            sol = bvp5c(@(z,y) base_ode(z,y,pk), @base_bc, solinit, opts);
        else
            sol = bvp4c(@(z,y) base_ode(z,y,pk), @base_bc, solinit, opts);
        end

        if p.verbose
            fprintf('mesh points = %d\n', numel(sol.x));
        end
    end
end

% ======================================================================= %
function dydz = base_ode(z,y,p)
    lambda = y(1,:);
    A = mobility_A(lambda,p.a);
    tau = 1 - z;
    dydz = zeros(size(y));
    dydz(1,:) = y(2,:);
    dydz(2,:) = (p.Gamma.*lambda.*tau.*A - 1 + lambda)./p.kappa;
end

% ======================================================================= %
function res = base_bc(ya,yb)
    res = [ya(2); yb(2)];
end

% ======================================================================= %
function y = base_initial_guess(z,p)
    % Initial guess based on the non-jammed algebraic branch where it exists.
    % Let r=1-lambda.  The diffusionless base equation has either r=0
    % (jammed branch) or Gamma*tau*(1-r)*(a+(1-a)r)=1.
    tau = max(1-z,0);
    r = zeros(size(z));
    for k = 1:numel(z)
        g0 = p.Gamma*tau(k)*p.a - 1;
        if g0 > 0
            lo = 0; hi = 1;
            for it = 1:60
                mid = 0.5*(lo+hi);
                gm = p.Gamma*tau(k)*(1-mid)*(p.a+(1-p.a)*mid) - 1;
                if gm > 0
                    lo = mid;
                else
                    hi = mid;
                end
            end
            r(k) = 0.5*(lo+hi);
        else
            r(k) = 0;
        end
    end
    lambda = 1 - r;

    if numel(z) > 1
        lambdap = gradient(lambda,z);
    else
        lambdap = 0*z;
    end
    y = [lambda; lambdap];
end

% ======================================================================= %
function kappas = continuation_sequence(kstart, kend, nsteps)
    if kend >= kstart
        kappas = kend;
        return;
    end
    if isempty(nsteps)
        decades = abs(log10(kstart/kend));
        nsteps = max(4, ceil(4*decades) + 1);
    end
    nsteps = max(2,nsteps);
    kappas = logspace(log10(kstart), log10(kend), nsteps);
end

% ======================================================================= %
function z = clustered_mesh(N,kappa)
    % Mesh for the nonlinear base BVP.  It clusters near z=1 where the
    % structural transition layer has thickness O(sqrt(kappa)).
    N = max(N,80);
    t = linspace(0,1,N);
    z1 = t;
    z2 = 1 - (1-t).^2;
    layer = min(0.35, max(8*sqrt(kappa), 0.02));
    z3 = 1 - layer*(1-t).^3;
    z3 = z3(z3>=0 & z3<=1);
    z = unique(sort([z1,z2,z3]));
    z(1) = 0;
    z(end) = 1;
end

% ======================================================================= %
function base = build_base_on_cheb_grid(p, sol)
    [z,D,D2] = cheb_on_01(p.N);
    Iop = integration_operator(D);

    vals = deval(sol,z);
    Lam = vals(1,:).';
    Lz  = vals(2,:).';
    A   = mobility_A(Lam,p.a);
    Ap  = mobility_A_prime(Lam,p.a);
    tau = 1 - z;
    Up = tau.*A;
    U  = integrate_from_zero(Iop, Up);
    Upp = -A + tau.*Ap.*Lz;
    Lzz = (p.Gamma.*Lam.*tau.*A - 1 + Lam)./p.kappa;
    Mdiag = 1./max(A,p.mobilityFloor);

    base = struct('N',p.N,'z',z,'D',D,'D2',D2,'Iop',Iop, ...
                  'Lambda',Lam,'Lambda_z',Lz,'Lambda_zz',Lzz, ...
                  'A',A,'A_lambda',Ap,'Mdiag',Mdiag, ...
                  'Up',Up,'U',U,'Upp',Upp,'sol',sol);

    if p.verbose
        fprintf('  base on Chebyshev grid: U_s=%.16g, Lambda=[%.6g, %.6g], min A=%.3e, max diagnostic M=%.3e\n', ...
            U(end), min(Lam), max(Lam), min(A), max(Mdiag));
    end
end

% ======================================================================= %
function lead = solve_leading_problem(p,base)
    if p.verbose
        fprintf('\nSolving leading O(1) long-wave problem...\n');
    end
    z = base.z; D = base.D; D2 = base.D2; Iop = base.Iop;
    Lam = base.Lambda; A = base.A; Ap = base.A_lambda;
    Up = base.Up; tau = 1-z;

    Lop = p.kappa*D2 - diag(1 + p.Gamma*Up + p.Gamma*Lam.*tau.*Ap);
    rhs = p.Gamma*Lam.*A;

    Mat = Lop;
    b = rhs;
    Mat(1,:) = D(1,:);       b(1) = 0;                  % ell0'(0)=0
    Mat(end,:) = D(end,:);   b(end) = -base.Lambda_zz(end); % ell0'(1)+Lambda''(1)=0

    ell0 = Mat \ b;
    u0p = A + tau.*Ap.*ell0;
    u0  = integrate_from_zero(Iop,u0p);
    psi0 = integrate_from_zero(Iop,u0);
    c = base.U(end) + psi0(end);

    lead = struct();
    lead.ell0 = ell0;
    lead.u0p = u0p;
    lead.u0 = u0;
    lead.psi0 = psi0;
    lead.c = c;

    if p.verbose
        fprintf('  c = %.16g, psi0(1)=%.16g\n', c, psi0(end));
    end
end

% ======================================================================= %
function corr = solve_correction_problem(p,base,lead,qFr)
    % Correct O(alpha) problem.  We solve the pressure-eliminated/vorticity
    % correction for the shear-stress coefficient T1.  This is essential for
    % recovering the Benjamin--Yih Newtonian threshold.
    z = base.z; D = base.D; D2 = base.D2; Iop = base.Iop;
    Lam = base.Lambda; Lz = base.Lambda_z;
    A = base.A; Ap = base.A_lambda; tau = 1-z;
    Up = base.Up; U = base.U; Upp = base.Upp;
    u0 = lead.u0; u0p = lead.u0p; psi0 = lead.psi0; ell0 = lead.ell0; c = lead.c;

    % T1'' = (FrV^2/S0)*[(U-c)*u0' - U''*psi0].
    Rvort = (U-c).*u0p - Upp.*psi0;
    MatT = D2;
    bT = (qFr/p.S0)*Rvort;
    MatT(end,:) = 0; MatT(end,end) = 1;          bT(end) = 0; % T1(1)=0
    MatT(end-1,:) = D(end,:);                    bT(end-1) = (qFr*(U(end)-c)*u0(end) + 1)/p.S0;
    T1 = MatT \ bT;
    T1p = D*T1;

    rlam = p.S0/p.T;
    if rlam <= 0
        error('Require rlambda=S0/T>0.');
    end

    Lop = p.kappa*D2 - diag(1 + p.Gamma*Up + p.Gamma*Lam.*tau.*Ap);
    rhs = ((U-c).*ell0 - psi0.*Lz)/rlam + p.Gamma*Lam.*A.*T1;

    MatM = Lop;
    bM = rhs;
    MatM(1,:) = D(1,:);       bM(1) = 0;       % m'(0)=0
    MatM(end,:) = D(end,:);   bM(end) = 0;     % m'(1)=0 with eta1=0 normalization
    m = MatM \ bM;

    vp = A.*T1 + tau.*Ap.*m;
    v = integrate_from_zero(Iop,vp);
    psi1 = integrate_from_zero(Iop,v);
    sigma2 = psi1(end);

    corr = struct();
    corr.qFr = qFr;
    corr.T1 = T1;
    corr.T1p = T1p;
    corr.m = m;
    corr.vp = vp;
    corr.v = v;
    corr.psi1 = psi1;
    corr.sigma2 = sigma2;
end

% ======================================================================= %
function A = mobility_A(lambda,a)
    % A(lambda)=1/M(lambda)=(1-lambda)*[a+(1-a)*(1-lambda)].
    r = max(1-lambda,0);
    A = r .* (a + (1-a).*r);
end

% ======================================================================= %
function Ap = mobility_A_prime(lambda,a)
    % dA/dlambda.  For lambda close to 1 this remains finite: Ap -> -a.
    r = max(1-lambda,0);
    Ap = -(a + 2*(1-a).*r);
end

% ======================================================================= %
function [z,D,D2] = cheb_on_01(N)
    % Chebyshev-Lobatto differentiation matrix on z in [0,1], ascending.
    if N < 2
        error('N must be at least 2.');
    end
    k = (0:N)';
    x = cos(pi*k/N);
    c = [2; ones(N-1,1); 2].*(-1).^k;
    X = repmat(x,1,N+1);
    dX = X - X';
    D_x = (c*(1./c)')./(dX + eye(N+1));
    D_x = D_x - diag(sum(D_x,2));

    z_desc = (x + 1)/2;
    D_z_desc = 2*D_x;
    P = flipud(eye(N+1));
    z = flipud(z_desc);
    D = P*D_z_desc*P;
    D2 = D*D;
end

% ======================================================================= %
function Iop = integration_operator(D)
    % Matrix approximation of F(z)=int_0^z f(s) ds using the same Chebyshev
    % differentiation matrix.  It solves D*F=f at collocation rows 2:end and
    % imposes F(0)=0 at row 1.
    n = size(D,1);
    B = D;
    B(1,:) = 0;
    B(1,1) = 1;
    Iop = zeros(n,n);
    for j = 1:n
        f = zeros(n,1);
        f(j) = 1;
        rhs = f;
        rhs(1) = 0;
        Iop(:,j) = B \ rhs;
    end
end

% ======================================================================= %
function F = integrate_from_zero(Iop,f)
    F = Iop*f(:);
end
