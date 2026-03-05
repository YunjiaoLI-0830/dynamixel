% function opt = optimize_depths_fig6(seg, kappa_des, l, epsFS, phi_o, g_bounds)
% %OPTIMIZE_DEPTHS_FIG6  Section III-B (Fig.6 Step2): optimize notch depths.
% %
% % Implements the paper's optimization described around Eq.(21):
% %   - Use Eq.(19) to enforce the desired curvature profile via compliance ratios
% %   - Use Eq.(20) (equal flexural rigidities) to reduce unknowns to g_o only
% %   - Choose the benchmark notch pair m such that max strain is highest, and
% %     solve g_o,m so that max strain equals epsFS (Eqs.(16)-(17))
% %   - Determine all other notch depths via alpha_j = alpha_m * (kappa_j/kappa_m)
% %
% % Inputs
% %   seg      : CPPRSegment (fields od_o,id_o,od_i,id_i,E and method calculate_section_properties)
% %   kappa_des: 1xn desired centerline curvatures (from Step1)
% %   l        : 1xn desired centerline bending lengths (Eq.14)
% %   epsFS    : allowable strain (epsilon_FS)
% %   phi_o    : 1xn notch orientations for outer tube (Eq.22)
% %   g_bounds : [gmin gmax] bounds for g_o (mm)
% %
% % Output struct opt:
% %   .g_o, .g_i (1xn)
% %   .I_o, .I_i, .gamma_o, .gamma_i (1xn)
% %   .d (1xn) distance between neutral axes
% %   .alpha (1xn) compliance constants (Eq.18)
% %   .m benchmark notch index used in final iteration
% %   .max_strain max strain across all notches at desired curvature
% %
% % Notes
% %   - This function uses the small-curvature approximation in Section III-B:
% %     kappa_o ~= kappa_i ~= kappa (centerline curvature).
% %
% 
% n = numel(kappa_des);
% assert(numel(phi_o)==n, 'phi_o must be 1xn');
% assert(numel(l)==n, 'l must be 1xn');
% 
% gmin = g_bounds(1);
% gmax = g_bounds(2);
% 
% % ---------- helpers ----------
%     function gi = gi_from_go(go)
%         gi = GeometryUtils.solve_gi_equal_EI(go, seg); % Eq.(20)
%     end
% 
%     function [Io, Ii, gam_o, gam_i, d, alpha] = props_from_go(go, ph)
%         gi = gi_from_go(go);
%         [Io, gam_o, Ii, gam_i] = seg.calculate_section_properties(go, gi, ph);
%         d = gam_o - gam_i;
%         alpha = d / (seg.E*Io + seg.E*Ii); % Eq.(18)
%     end
% 
%     function epsm = max_strain_pair(go, kappa, ph)
%         gi = gi_from_go(go);
%         epsm = GeometryUtils.max_strain_notch(kappa, go, gi, ph, seg); % Eq.(16)-(17)
%     end
% 
%     function go = solve_go_for_eps(kappa, ph)
%         % Solve go such that max strain in this pair equals epsFS
%         f = @(go) max_strain_pair(go, kappa, ph) - epsFS;
%         % Reasonable initial guess (paper's initial guess uses Eq.(16) edge a)
%         go0 = min(max(0.25*(seg.od_o/2), gmin), gmax);
%         go = GeometryUtils.robust_fzero_scalar(f, go0, gmin, gmax);
%     end
% 
%     function go = solve_go_for_alpha(alpha_target, ph)
%         % Solve go such that alpha(go) matches alpha_target (with gi from Eq.(20))
%         f = @(go) (props_from_go(go, ph)); %#ok<NASGU>
%         % MATLAB doesn't allow multiple outputs in anonymous; implement inline:
%         g = @(go) alpha_from_go(go, ph) - alpha_target;
%         go0 = min(max(0.25*(seg.od_o/2), gmin), gmax);
%         go = GeometryUtils.robust_fzero_scalar(g, go0, gmin, gmax);
%     end
% 
%     function a = alpha_from_go(go, ph)
%         [~, ~, ~, ~, ~, a] = props_from_go(go, ph);
%     end
% 
% % ---------- initialization ----------
% % Choose initial benchmark m as highest |kappa|
% [~, m] = max(abs(kappa_des));
% if abs(kappa_des(m)) < 1e-12
%     % all curvatures ~0 => shallow notches (max stiffness)
%     opt = pack_constant_depths(gmin*ones(1,n));
%     return;
% end
% 
% % initialize go pattern proportional to |kappa| (paper suggestion)
% go = zeros(1,n);
% go(m) = min(max(0.25*(seg.od_o/2), gmin), gmax);
% for j = 1:n
%     go(j) = go(m) * abs(kappa_des(j)/kappa_des(m));
%     go(j) = min(max(go(j), gmin), gmax);
% end
% 
% % ---------- iterate 1-3 times to ensure benchmark selection is correct ----------
% for outer_iter = 1:3
%     % (1) solve benchmark depth so that its max strain hits epsFS
%     go(m) = solve_go_for_eps(kappa_des(m), phi_o(m));
% 
%     % (2) compute benchmark alpha_m
%     gi_m = gi_from_go(go(m));
%     [alpha_m, ~, ~, ~, ~, ~] = GeometryUtils.compliance_alpha(go(m), gi_m, phi_o(m), seg);
% 
%     % (3) enforce Eq.(19): alpha_j = alpha_m * (kappa_j/kappa_m)
%     for j = 1:n
%         if j == m, continue; end
%         if abs(kappa_des(j)) < 1e-12
%             go(j) = gmin; % near-straight segment -> max stiffness
%             continue;
%         end
%         alpha_target = alpha_m * (kappa_des(j) / kappa_des(m));
%         % alpha should have same sign as kappa (since tau is scalar)
%         % our alpha(...) already encodes sign through d (via phi_o)
%         go(j) = solve_go_for_alpha(alpha_target, phi_o(j));
%     end
% 
%     % (4) evaluate max strain across all notches; update benchmark if needed
%     eps_all = zeros(1,n);
%     for j = 1:n
%         eps_all(j) = max_strain_pair(go(j), kappa_des(j), phi_o(j));
%     end
%     [eps_max, m_new] = max(eps_all);
% 
%     % if benchmark already the max strain (or ties), stop
%     if m_new == m
%         break;
%     end
%     m = m_new;
% end
% 
% % ---------- finalize properties ----------
% g_o = go;
% g_i = zeros(1,n);
% I_o = zeros(1,n); I_i = zeros(1,n);
% gamma_o = zeros(1,n); gamma_i = zeros(1,n);
% d = zeros(1,n); alpha = zeros(1,n);
% eps_all = zeros(1,n);
% 
% for j = 1:n
%     g_i(j) = gi_from_go(g_o(j));
%     [I_o(j), gamma_o(j), I_i(j), gamma_i(j)] = seg.calculate_section_properties(g_o(j), g_i(j), phi_o(j));
%     d(j) = gamma_o(j) - gamma_i(j);
%     alpha(j) = d(j) / (seg.E*I_o(j) + seg.E*I_i(j));
%     eps_all(j) = GeometryUtils.max_strain_notch(kappa_des(j), g_o(j), g_i(j), phi_o(j), seg);
% end
% 
% opt = struct();
% opt.g_o = g_o;
% opt.g_i = g_i;
% opt.I_o = I_o;
% opt.I_i = I_i;
% opt.gamma_o = gamma_o;
% opt.gamma_i = gamma_i;
% opt.d = d;
% opt.alpha = alpha;
% opt.m = m;
% opt.eps_per_notch = eps_all;
% opt.max_strain = max(eps_all);
% opt.epsFS = epsFS;
% 
% end
% 
% function opt = pack_constant_depths(go)
% opt = struct();
% opt.g_o = go;
% opt.g_i = go;
% opt.I_o = nan(size(go));
% opt.I_i = nan(size(go));
% opt.gamma_o = nan(size(go));
% opt.gamma_i = nan(size(go));
% opt.d = nan(size(go));
% opt.alpha = nan(size(go));
% opt.m = 1;
% opt.eps_per_notch = zeros(size(go));
% opt.max_strain = 0;
% opt.epsFS = 0;
% end


%-------------------------------new version--------------


function opt = optimize_depths_fig6(kappa, seg, epsFS, g_bounds, phi_o)
%OPTIMIZE_DEPTHS_FIG6  Paper Sec.III-B notch depth optimization (Fig.6 step 2).
%
% Inputs
%   kappa   1xn desired centerline curvatures κ_j
%   seg     CPPRSegment (has od_o,id_o,od_i,id_i,E and calculate_section_properties)
%   epsFS   strain limit ε_FS
%   g_bounds [gmin gmax] bounds for go (mm)
%   phi_o   1xn notch orientations (can pass zeros or use sign rule later)
%
% Output (struct)
%   opt.go, opt.gi, opt.alpha, opt.gamma_o, opt.gamma_i, opt.eps_max, opt.m

n = numel(kappa);
gmin = g_bounds(1); gmax = g_bounds(2);

% ---- choose initial benchmark (highest |kappa|) ----
[~, m] = max(abs(kappa));
if abs(kappa(m)) < 1e-12
    error('All desired curvatures are ~0; depth optimization is ill-posed.');
end

% ---- outer iteration: update benchmark if needed (paper says 1-2 rounds) ----
for outer = 1:3
    % Solve go_m such that max strain hits epsFS
    go_m = fzero(@(go) strain_margin_for_benchmark(go, m), clamp_mid(gmin,gmax));
    go_m = min(max(go_m, gmin), gmax);

    % Build full depth pattern from this benchmark go_m
    [go, gi, alpha, gam_o, gam_i, eps_max] = pattern_from_benchmark(go_m, m);

    % verify true highest-strained notch
    [~, m_new] = max(eps_max);
    if m_new == m
        % converged benchmark
        opt = struct();
        opt.go = go; opt.gi = gi;
        opt.alpha = alpha;
        opt.gamma_o = gam_o; opt.gamma_i = gam_i;
        opt.eps_max = eps_max;
        opt.m = m;
        opt.go_m = go_m;
        opt.outer_iters = outer;
        return;
    end

    % update benchmark and repeat
    m = m_new;
end

warning('Benchmark update did not converge in 3 iterations; returning last result.');
opt = struct();
opt.go = go; opt.gi = gi;
opt.alpha = alpha;
opt.gamma_o = gam_o; opt.gamma_i = gam_i;
opt.eps_max = eps_max;
opt.m = m;
opt.go_m = go_m;
opt.outer_iters = outer;

% ================= nested helpers =================

    function val = strain_margin_for_benchmark(go_trial, m_idx)
        % For a candidate benchmark go_m, compute full pattern, then margin to epsFS:
        % val = max(eps) - epsFS  (want zero)
        go_trial = min(max(go_trial, gmin), gmax);
        [~, ~, ~, ~, ~, eps_max_trial] = pattern_from_benchmark(go_trial, m_idx);
        val = max(eps_max_trial) - epsFS;
    end

    function [go, gi, alpha, gam_o, gam_i, eps_max] = pattern_from_benchmark(go_m_local, m_idx)
        go = zeros(1,n);
        gi = zeros(1,n);
        alpha = zeros(1,n);
        gam_o = zeros(1,n);
        gam_i = zeros(1,n);
        eps_max = zeros(1,n);

        % First compute benchmark properties
        gi_m = GeometryUtils.solve_gi_equal_EI(go_m_local, seg);
        [alpha_m, ~, ~, ~, gam_o_m, gam_i_m] = GeometryUtils.compliance_alpha(go_m_local, gi_m, phi_o(m_idx), seg);

        % Desired alpha ratios from Eq.(19): alpha_j = alpha_m * (kappa_j/kappa_m)
        kappa_m = kappa(m_idx);

        % Solve each notch independently for go_j (1D root) so that alpha_j matches target
        for j = 1:n
            target_alpha_j = alpha_m * (kappa(j) / kappa_m);

            if abs(kappa(j)) < 1e-12
                % near straight: choose minimal cut (stiffest)
                go(j) = gmin;
                gi(j) = GeometryUtils.solve_gi_equal_EI(go(j), seg);
            else
                % Use paper's suggested initial guess: go_j ~ go_m * |kappa_j/kappa_m|
                go0 = go_m_local * abs(kappa(j) / kappa_m);
                go0 = min(max(go0, gmin), gmax);

                % Root function: alpha(go) - target_alpha_j = 0
                f = @(goj) alpha_of_go(goj, j) - target_alpha_j;

                % Robust solve in [gmin, gmax]
                go(j) = GeometryUtils.robust_fzero_scalar(f, go0, gmin, gmax);
                gi(j) = GeometryUtils.solve_gi_equal_EI(go(j), seg);
            end

            % record alpha + gammas + strain
            [alpha(j), ~, ~, ~, gam_o(j), gam_i(j)] = GeometryUtils.compliance_alpha(go(j), gi(j), phi_o(j), seg);
            eps_max(j) = GeometryUtils.max_strain_notch(kappa(j), go(j), gi(j), phi_o(j), seg);
        end

        % overwrite benchmark with the exact local values (to avoid tiny numerical drift)
        go(m_idx) = go_m_local;
        gi(m_idx) = gi_m;
        alpha(m_idx) = alpha_m;
        gam_o(m_idx) = gam_o_m;
        gam_i(m_idx) = gam_i_m;
        eps_max(m_idx) = GeometryUtils.max_strain_notch(kappa(m_idx), go_m_local, gi_m, phi_o(m_idx), seg);
    end

    function a = alpha_of_go(goj, j)
        goj = min(max(goj, gmin), gmax);
        gij = GeometryUtils.solve_gi_equal_EI(goj, seg);
        [a, ~] = GeometryUtils.compliance_alpha(goj, gij, phi_o(j), seg);
    end

end

function x = clamp_mid(lo, hi)
x = 0.5*(lo+hi);
end
