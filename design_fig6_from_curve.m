function out = design_fig6_from_curve(curve, params)
%DESIGN_FIG6_FROM_CURVE  Fig.6 design flow (Section III) implementation.
%
% This ties together:
%   Step 1 (Section III-A): sample desired curve -> {kappa_j, l_j, c}
%   Step 2 (Section III-B): optimize notch depths (go,gi) under strain limit
%   Step 3 (Section III-C): compute remaining parameters (phi_o, h) using Eq.(10),(22)
%
% Inputs
%   curve  struct describing the desired curve (see sample_curve_for_fig6)
%   params struct fields (typical):
%     .Lcurve   total desired curve length (mm)
%     .beta     straight fraction (Eq.14)
%     .n        number of notches
%     .epsFS    strain limit with factor of safety (epsilon_FS)
%     .segment_params  CPPRSegment object (od_o,id_o,od_i,id_i,E)
%     (optional) .g_bounds [gmin gmax] bounds for g_o (mm)
%     (optional) .method   passed to sample_curve_for_fig6 ('sequential' or 'integral')
%
% Output out:
%   .notches : struct with fields h,g_o,g_i,phi_o,c
%   .sample  : output from sample_curve_for_fig6
%   .opt     : output from optimize_depths_fig6
%   .debug   : optional forward-model check using CPPRRobot
%

arguments
    curve struct
    params struct
end

% --- robust handling for qmax (Table II uses qmax_mm) ---
if ~isfield(params, 'qmax') || isempty(params.qmax)
    if isfield(params, 'qmax_mm'), params.qmax = params.qmax_mm; end
    if isfield(params, 'qMax'),    params.qmax = params.qMax;    end
end
if ~isfield(params,'qmax') || isempty(params.qmax)
    % Default: design at qmax = 0 (still allows computing notch geometry)
    params.qmax = 0;
end


seg = params.segment_params;

% defaults
if ~isfield(params,'g_bounds')
    gmin = 0.2;
    gmax = max(gmin+0.1, seg.od_o/2 - 0.2);
    params.g_bounds = [gmin gmax];
end
if ~isfield(params,'method')
    params.method = "sequential";
end

% ---------------- Step 1 sample ----------------
sample = sample_curve_for_fig6(curve, params.Lcurve, params.n, params.beta, params.method);
kappa = sample.kappa(:)';   % 1xn
l = sample.l(:)';           % 1xn (centerline bend lengths)
c = sample.c;               % scalar (rigid length per module)

%-----------------Step A Strain Modeling-----------------------------
% GeometryUtils.max_strain_notch
% uncertain!!!!!!!!!

% ---------------- Step 3a (III-C): phi_o via Eq.(22) ----------------
% Eq.(22): phi_o,j = 0 if kappa_j < 0, pi if kappa_j >= 0
phi_o = zeros(1, params.n);
phi_o(kappa >= 0) = pi;

% ---------------- Step 2 B optimize depths ----------------
% opt = optimize_depths_fig6(seg, kappa, l, params.epsFS, phi_o, params.g_bounds);
opt = optimize_depths_fig6(kappa, seg, params.epsFS, params.g_bounds, phi_o);

% ---------------- Step 3C notch lengths h via Eq.(10) ----------------
% Eq.(10): l_j = h_j / (1 - gamma_o,j * kappa_j)  ->  h_j = l_j * (1 - gamma_o,j * kappa_j)
h = zeros(1, params.n);
for j = 1:params.n
    factor = 1 - opt.gamma_o(j) * kappa(j);

    % numerical safety (should rarely trigger if design is valid)
    if factor <= 0
        error('Invalid geometry: 1 - gamma_o * kappa <= 0 at notch %d', j);
    end

    % Eq.(III-C): l_j = h_j * (1 - gamma_o,j * kappa_j)
    % => h_j = l_j * (1 - gamma_o,j * kappa_j)
    h(j) = l(j) * factor;
end

% paper-stated bounds:  l_j <= h_j < 2*l_j
h = min(max(h, l), 2*l - 1e-9);


% paper implies: l_j <= h_j < 2 l_j (self-interference at kappa = (OD/2)^-1)
h = max(h, l);
h = min(h, 2*l - 1e-9);

% ---------------- pack for CPPRRobot ----------------
notches = struct();
notches.h = h;
notches.g_o = opt.g_o;
notches.g_i = opt.g_i;
notches.phi_o = phi_o;
notches.c = c;

out = struct();
out.notches = notches;
out.sample = sample;
out.opt = opt;

% ---------------- optional validation using forward model ----------------
try
    robot = CPPRRobot(seg, notches);
    % Use Eq.(8)/(7) solver in your CPPRRobot (displacement-driven)
    % Here we compute the required q for the desired centerline curvature profile:
    %   q = sum_j d_j * h_j * kappa_o,j (Eq.8), but kappa_o,j differs from kappa_j.
    % For a quick check, we actuate the model using q that matches centerline curvature approx:
    q_est = sum(opt.d(:) .* h(:) .* kappa(:));
    kappa_o = robot.solve_shape(q_est);
    pts = robot.get_backbone_shape(kappa_o);

    out.debug.q_est = q_est;
    out.debug.kappa_o = kappa_o(:)';
    out.debug.tip = pts(:,end);
    out.debug.max_strain = opt.max_strain;
catch ME
    out.debug.error = ME.message;
end

end