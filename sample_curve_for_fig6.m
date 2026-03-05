function sample = sample_curve_for_fig6(curve, Lcurve, n, beta, method)
%SAMPLE_CURVE_FOR_FIG6  Fig.6 curve discretization (Sec. III) producing desired
%curvatures kappa_j for each notched segment.
%
% This file implements TWO ways to obtain kappa_j:
%   (A) 'sequential' (DEFAULT): matches the paper text:
%       starting at the base, solve ONE scalar nonlinear equation per step to
%       pick kappa_j that minimizes the normal error at the next sampled point.
%       (Uses transforms equivalent to Eqs. (11)–(13) for a constant-curvature
%       notch followed by a rigid segment.)
%   (B) 'integral' : integrates a known curvature function over the bend interval.
%
% Inputs
%   curve  struct:
%     - type='segments', curve.segments: array with fields .L (mm), .theta (rad)
%     - type='points',   curve.p: Nx2 or Nx3 points (mm)
%   Lcurve total desired arc-length (mm)
%   n      number of modules (notched-rigid pairs)
%   beta   straight fraction within each module (Eq.14)
%   method (optional) 'sequential' or 'integral'
%
% Outputs
%   sample.c       straight length per module (Eq.14): c = (L/n)*beta
%   sample.l       bend length per module     (Eq.14): l = (L/n)*(1-beta)
%   sample.s0,s1   module start/end arc-lengths
%   sample.p_des   3x(n+1) desired points at s = [0,s1]
%   sample.t_des   3x(n+1) desired tangents at s = [0,s1]
%   sample.n_des   3x(n+1) desired normals  at s = [0,s1] (planar)
%   sample.kappa   1xn desired curvature per notched segment
%   sample.theta   1xn desired turning of notched segment (theta_j = kappa_j*l)
%
% Notes
%  - Coordinate convention follows your CPPRRobot.get_backbone_shape:
%    the backbone grows along +z, bending in the x-z plane (y=0).
%

arguments
    curve struct
    Lcurve (1,1) double {mustBePositive}
    n (1,1) double {mustBeInteger,mustBePositive}
    beta (1,1) double {mustBeGreaterThanOrEqual(beta,0),mustBeLessThanOrEqual(beta,1)}
    method (1,1) string = "sequential"
end

% --- Eq.(14): module lengths ---
Lmod = Lcurve / n;
c = Lmod * beta;
l = Lmod * (1 - beta);

% --- module arc-length edges ---
s0 = (0:n-1) * Lmod;
s1 = (1:n)   * Lmod;

% In the paper, p_j is the end of the j-th notched-rigid pair -> sampled at s1(j)
s_samples = [0, s1];

% --- build desired curve evaluator: position + tangent at arc-length s ---
eval_des = build_desired_curve_evaluator(curve, Lcurve);

p_des = zeros(3, n+1);
t_des = zeros(3, n+1);
n_des = zeros(3, n+1);

% eq.15
for k=1:(n+1)
    [p_des(:,k), t_des(:,k)] = eval_des(s_samples(k));
    % planar normal: rotate tangent by +90deg about y (x-z plane)
    tx = t_des(1,k); tz = t_des(3,k);
    nn = [tz; 0; -tx];
    if norm(nn) < 1e-12
        nn = [1;0;0];
    end
    n_des(:,k) = nn / norm(nn);
end

theta = zeros(1,n);
kappa = zeros(1,n);

switch lower(method)
    case "integral"
        % Old-style: integrate known kappa(s) over the bend interval.
        % (Kept for debugging; paper-matching default is 'sequential'.)
        % Build piecewise-constant curvature function for 'segments', or numeric curvature for 'points'.
        [s_u, kappa_u] = build_curvature_function(curve, Lcurve, n);
        for j=1:n
            % bend interval within module: last l portion -> [s1-l, s1]
            a = s1(j) - l;
            b = s1(j);
            theta(j) = trapz_on_interval(s_u, kappa_u, a, b);
        end
        kappa = theta ./ max(l, eps);

    case "sequential"
        % Paper-matching: one nonlinear equation per step (normal error at next point).
        % Maintain current pose T in SE(2) embedded in 4x4 homogeneous matrix (x,z plane).
        T = eye(4);
        % We use the desired normal at s1(j) (end of pair) as the error direction.
        for j=1:n
            p_target = p_des(:, j+1);
            n_target = n_des(:, j+1);

            % Initial guess: local desired curvature over the *bend* length.
            % Approximate using change in desired tangent direction over the module / l.
            t_prev = t_des(:, j);
            t_next = t_des(:, j+1);
            dtheta_mod = atan2( t_prev(1)*t_next(3) - t_prev(3)*t_next(1), ...
                                t_prev(1)*t_next(1) + t_prev(3)*t_next(3) );
            k0 = dtheta_mod / max(l, eps);

            f = @(kj) normal_error_after_pair(T, kj, l, c, p_target, n_target);
            % Robust scalar solve: fzero around k0
            kappa(j) = robust_fzero(f, k0);

            theta(j) = kappa(j) * l;

            % Update pose T using the same constant-curvature transform form as CPPRRobot (bend then rigid)
            T = T * T_bend(kappa(j), l) * T_rigid(c);
        end

    otherwise
        error("Unknown method '%s'. Use 'sequential' or 'integral'.", method);
end

sample = struct();
sample.c = c;
sample.l = repmat(l, 1, n);
sample.s0 = s0;
sample.s1 = s1;
sample.p_des = p_des;
sample.t_des = t_des;
sample.n_des = n_des;
sample.kappa = kappa;
sample.theta = theta;
sample.method = method;

end

% ======================= helpers =======================

function val = normal_error_after_pair(Tprev, kj, l, c, p_target, n_target)
% Compute normal component of error at end of (bend+rigid) given current pose.
Tnext = Tprev * T_bend(kj, l) * T_rigid(c);
p_model = Tnext(1:3,4);
e = p_model - p_target;
val = dot(n_target, e);
end

function T = T_bend(k, L)
% Constant-curvature bend transform in x-z plane (matching CPPRRobot convention).
% If k ~ 0 -> straight along +z
if abs(k) < 1e-9
    T = [1 0 0 0;
         0 1 0 0;
         0 0 1 L;
         0 0 0 1];
else
    theta = k*L;
    s = sin(theta); c = cos(theta); inv_k = 1/k;
    T = [ c 0  s  inv_k*(1-c);
          0 1  0          0;
         -s 0  c    inv_k*s;
          0 0  0          1];
end
end

function T = T_rigid(L)
% Rigid straight segment along +z
T = [1 0 0 0;
     0 1 0 0;
     0 0 1 L;
     0 0 0 1];
end

function root = robust_fzero(f, k0)
% Try fzero with expanding brackets around k0; fallback to local minimization.
opts = optimset('Display','off');
try
    % First try: direct
    root = fzero(f, k0, opts);
    return;
catch
end
% Bracket search
span = 0.5;
for it=1:12
    a = k0 - span;
    b = k0 + span;
    try
        fa = f(a); fb = f(b);
        if sign(fa) ~= sign(fb)
            root = fzero(f, [a b], opts);
            return;
        end
    catch
    end
    span = span * 2;
end
% Fallback: minimize |f|
g = @(k) abs(f(k));
root = fminsearch(g, k0, opts);
end

function eval_des = build_desired_curve_evaluator(curve, Lcurve)
% Returns function handle eval_des(s)-> [p,t] (3x1 each), for s in [0,Lcurve].
switch lower(string(curve.type))
    case "segments"
        segs = curve.segments;
        if isempty(segs), error("curve.segments is empty."); end

        % Build segment list with lengths scaled to match Lcurve, and constant curvature per segment.
%         L_list = arrayfun(@(s) s.L, segs);
%         th_list = arrayfun(@(s) s.theta, segs);
          L_list  = arrayfun(@(s) s.L, segs);    L_list  = L_list(:).';   % 强制 1×m
          th_list = arrayfun(@(s) s.theta, segs); th_list = th_list(:).'; % 强制 1×m


        L_in = sum(L_list);
        if abs(L_in - Lcurve) > 1e-9
            scale = Lcurve / L_in;
            L_list = L_list * scale;
            th_list = th_list * scale; % keep curvature shape (k=theta/L) unchanged
        end
        k_list = th_list ./ max(L_list, eps);
        s_edges = [0, cumsum(L_list)];

        % Precompute pose at each segment boundary for fast eval
        T_edges = cell(1, numel(L_list)+1);
        T_edges{1} = eye(4);
        for i=1:numel(L_list)
            T_edges{i+1} = T_edges{i} * T_bend(k_list(i), L_list(i)); % whole segment is "bend"
        end

        eval_des = @(s) eval_piecewise_const_curve(s, s_edges, k_list, T_edges);

    case "points"
        P = curve.p;
        if size(P,2) > 2, P = P(:,[1 3]); end
        ds = sqrt(sum(diff(P,1,1).^2,2));
        s = [0; cumsum(ds)];
        if s(end) < 1e-9, error("Curve points have near-zero length."); end

        % scale s to Lcurve (so caller's s matches Table II length)
        scale = Lcurve / s(end);
        s = s * scale;
        P = P * scale; % uniform scaling in x-z

        % Build interpolants for position and tangent (numerical derivative)
        eval_des = @(su) eval_points_curve(su, s, P);

    otherwise
        error("Unsupported curve.type '%s'.", curve.type);
end
end

function [p,t] = eval_piecewise_const_curve(s, s_edges, k_list, T_edges)
s = min(max(s, s_edges(1)), s_edges(end));
% find segment index i s.t. s in [s_edges(i), s_edges(i+1)]
i = find(s_edges <= s, 1, 'last');
i = min(i, numel(k_list));
s0 = s_edges(i);
ds = s - s0;
T = T_edges{i} * T_bend(k_list(i), ds);
p = T(1:3,4);
% tangent is T * e_z
ez = [0;0;1;0];
tv = T * ez;
t = tv(1:3);
t = t / max(norm(t), 1e-12);
end

function [p,t] = eval_points_curve(su, s, P)
% Linear interp for p, central diff for t
x = interp1(s, P(:,1), su, 'linear','extrap');
z = interp1(s, P(:,2), su, 'linear','extrap');
p = [x; 0; z];

% tangent from finite differences on interpolated curve
du = 1e-3 * (s(end)-s(1)); % small step
s1 = min(s(end), su+du);
s0 = max(s(1), su-du);
x1 = interp1(s, P(:,1), s1, 'linear','extrap');
z1 = interp1(s, P(:,2), s1, 'linear','extrap');
x0 = interp1(s, P(:,1), s0, 'linear','extrap');
z0 = interp1(s, P(:,2), s0, 'linear','extrap');
t = [x1-x0; 0; z1-z0];
t = t / max(norm(t), 1e-12);
end

function [s_u, kappa_u] = build_curvature_function(curve, Lcurve, n)
% For 'integral' method.
switch lower(string(curve.type))
    case "segments"
        segs = curve.segments;
        L_list = arrayfun(@(s) s.L, segs);
        th_list = arrayfun(@(s) s.theta, segs);
        L_in = sum(L_list);
        if abs(L_in - Lcurve) > 1e-9
            scale = Lcurve / L_in;
            L_list = L_list * scale;
            th_list = th_list * scale;
        end
        k_list = th_list ./ max(L_list, eps);
        s_edges = [0, cumsum(L_list)];
        % represent as dense samples for generic trapz
        m = max(2000, 50*n);
        s_u = linspace(0, Lcurve, m)';
        kappa_u = zeros(size(s_u));
        for i=1:numel(k_list)
            mask = (s_u >= s_edges(i)) & (s_u <= s_edges(i+1));
            kappa_u(mask) = k_list(i);
        end
    case "points"
        P = curve.p;
        if size(P,2) > 2, P = P(:,[1 3]); end
        ds = sqrt(sum(diff(P,1,1).^2,2));
        s = [0; cumsum(ds)];
        scale = Lcurve / s(end);
        P = P * scale;
        s = s * scale;
        m = max(2000, 50*n);
        s_u = linspace(0, Lcurve, m)';
        x = interp1(s, P(:,1), s_u, 'linear');
        z = interp1(s, P(:,2), s_u, 'linear');
        ds_u = s_u(2)-s_u(1);
        x1 = gradient(x, ds_u); z1 = gradient(z, ds_u);
        x2 = gradient(x1, ds_u); z2 = gradient(z1, ds_u);
        denom = (x1.^2 + z1.^2).^(3/2);
        kappa_u = zeros(size(denom));
        mask = denom > 1e-12;
        kappa_u(mask) = (x1(mask).*z2(mask) - z1(mask).*x2(mask)) ./ denom(mask);
    otherwise
        error("Unsupported curve.type for curvature function.");
end
end

function th = trapz_on_interval(s_u, kappa_u, a, b)
a = max(a, s_u(1));
b = min(b, s_u(end));
mask = (s_u>=a) & (s_u<=b);
if nnz(mask) < 2
    th = 0;
else
    th = trapz(s_u(mask), kappa_u(mask));
end
end
