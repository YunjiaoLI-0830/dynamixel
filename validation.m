%% CPPR Validation & Design Script
% Based on Oliver-Butler 2022 - Concentric Push-Pull Robots Planar Modeling and Design
clear; clc; close all;

%% ================= 1. Parameters & Initialization =================
% System Parameters
n = 10;                 % Number of notches
E = 75e3;               % Young's Modulus (N/mm^2) - Nitinol
OD_o = 4.02; ID_o = 3.70; % Outer Tube Geometry (Assume thin wall, specific IDs needed)
OD_i = 3.50; ID_i = 3.20; % Inner Tube Geometry
% Note: Paper uses wall thicknesses, here approximated by ID = OD - 2*t
% t_o = 0.16 -> ID_o = 3.70; t_i = 0.15 -> ID_i = 3.20

% Design Constraints
epsilon_max = 0.02;     % Material strain limit (2%)
FS = 1.5;               % Factor of Safety
epsilon_FS = epsilon_max / FS; 

% Shape Definition
beta = 0.5;             % Paper suggests beta approx 0.5 (Eq. 14 discussion)
L_total = 100;          % Total length of robot (mm)
l_segment = (L_total / n) * (1 - beta); % Length of curved section
c = (L_total / n) * beta;               % Length of rigid section (spacing)

% Target Curvatures (Sampled from desired shape)
kappa_desired = [0.01, 0.02, 0.03, 0.025, 0.015, -0.015, -0.025, -0.03, -0.02, -0.01];

%% ================= 2. Benchmark Notch Design (Section III-B) =================
fprintf('--- Starting Design Optimization ---\n');

% Find the notch with maximum curvature magnitude
[max_kappa, m_idx] = max(abs(kappa_desired));
fprintf('Benchmark Notch Index: %d, Max Curvature: %.4f mm^-1\n', m_idx, max_kappa);

% --- Helper Functions for Geometry ---
% We need functions that calculate I (Inertia) and gamma (Neutral Axis) 
% based on cut depth g. See 'geometry_properties' at end of script.

% Optimization: Find g_o for benchmark pair that hits strain limit
% Objective: max(strain_edge_a, strain_edge_b) - epsilon_FS = 0
options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);

% Initial guess for cut depth (e.g., radius)
g_guess = OD_o / 2; 

% Solve for Outer Tube Depth (g_o) that satisfies strain limit
solve_strain = @(g) calculate_max_strain(g, max_kappa, OD_o, ID_o) - epsilon_FS;
g_o_bench = fsolve(solve_strain, g_guess, options);

% Solve for matching Inner Tube Depth (g_i) (Eq. 20: EI_o = EI_i)
% Since E is same, we match I_o = I_i
[I_o_bench, ~, ~] = geometry_properties(OD_o, ID_o, g_o_bench);
solve_stiffness = @(g) geometry_properties(OD_i, ID_i, g) - I_o_bench; % Match Inertia
g_i_bench = fsolve(solve_stiffness, OD_i/2, options);

fprintf('Benchmark Found: g_o = %.3f mm, g_i = %.3f mm\n', g_o_bench, g_i_bench);

%% ================= 3. Pattern Generation (Eq. 21) =================
% Calculate Compliance Constant alpha for benchmark (Eq. 18)
% alpha = d / (EI_o + EI_i)
[~, gamma_o_bench, ~] = geometry_properties(OD_o, ID_o, g_o_bench);
[~, gamma_i_bench, ~] = geometry_properties(OD_i, ID_i, g_i_bench);
% Note: In design phase, we assume aligned neutral axes for d calc
d_bench = abs(gamma_i_bench) + abs(gamma_o_bench); 
EI_sum_bench = E * (I_o_bench + I_o_bench); % I_o = I_i enforced
alpha_bench = d_bench / EI_sum_bench;

g_o_pattern = zeros(1, n);
g_i_pattern = zeros(1, n);
phi_pattern = zeros(1, n);

for j = 1:n
    k_j = kappa_desired(j);
    
    % 1. Determine Orientation (Eq. 22)
    if k_j >= 0
        phi_pattern(j) = pi; % Outer tube notch at pi (pulling inner)
    else
        phi_pattern(j) = 0;  % Outer tube notch at 0
    end
    
    if abs(k_j) < 1e-6
        % Handle straight sections (minimal cut or zero curvature)
        g_o_pattern(j) = 0; 
        g_i_pattern(j) = 0;
        continue;
    end

    % 2. Solve for g_o using Compliance Scaling (Eq. 21)
    % Constraint: alpha_j = (alpha_m / k_m) * k_j
    target_alpha = (alpha_bench / abs(max_kappa)) * abs(k_j);
    
    % Function to solve: calc_alpha(g) - target_alpha = 0
    % We must solve g_o and g_i simultaneously or nested to maintain EI_o = EI_i
    solve_shape = @(g_o) solve_alpha_discrepancy(g_o, target_alpha, E, OD_o, ID_o, OD_i, ID_i) ;
    g_o_pattern(j) = fsolve(solve_shape, g_o_bench * abs(k_j/max_kappa), options);
    
    % 3. Calculate corresponding g_i (Stiffness Matching)
    [I_o_val, ~, ~] = geometry_properties(OD_o, ID_o, g_o_pattern(j));
    solve_stiffness_j = @(g) geometry_properties(OD_i, ID_i, g) - I_o_val;
    g_i_pattern(j) = fsolve(solve_stiffness_j, OD_i/2, options);
end

%% ================= 4. Calculate Notch Lengths (Eq. 10) =================
h_pattern = zeros(1, n);
gamma_o_vec = zeros(1, n); % Signed gamma
gamma_i_vec = zeros(1, n); % Signed gamma

for j = 1:n
    [~, gam_o_mag, ~] = geometry_properties(OD_o, ID_o, g_o_pattern(j));
    [~, gam_i_mag, ~] = geometry_properties(OD_i, ID_i, g_i_pattern(j));
    
    % Apply orientation (Eq. 1)
    % gamma_o is usually negative of gamma_i direction relative to center
    gamma_o_vec(j) = gam_o_mag * cos(phi_pattern(j));
    gamma_i_vec(j) = gam_i_mag * cos(phi_pattern(j) - pi);
    
    % Eq. 10: l_j = h_j / (1 - gamma_o * kappa_j)
    % => h_j = l_segment * (1 - gamma_o * kappa)
    % Note: kappa here is centerline curvature
    % Paper uses Eq 10 to find l_j from h, or h from l_j. 
    % We fixed l_segment (desired arc length), so we solve for h.
    h_pattern(j) = l_segment * (1 - gamma_o_vec(j) * kappa_desired(j));
end

%% ================= 5. Forward Kinematics Verification =================
fprintf('--- Running Forward Kinematics Verification ---\n');

% Calculate Geometric Properties for all notches
I_o_vec = zeros(1,n); I_i_vec = zeros(1,n);
d_vec = zeros(1,n);

for j = 1:n
    [I_o_vec(j), ~, ~] = geometry_properties(OD_o, ID_o, g_o_pattern(j));
    [I_i_vec(j), ~, ~] = geometry_properties(OD_i, ID_i, g_i_pattern(j));
    
    % Eq. 3: d = gamma_i - gamma_o
    d_vec(j) = gamma_i_vec(j) - gamma_o_vec(j);
end

% Calculate Required Actuation q (Eq. 8)
% q = sum( d_j * h_j * kappa_o,j )
% Note: kappa_desired is Centerline curvature. We need Outer Backbone curvature for Eq 8.
% Eq. 9 Inverse: kappa_o = kappa / (1 - gamma_o * kappa)
kappa_o_desired = kappa_desired ./ (1 - gamma_o_vec .* kappa_desired);
q_required = sum(d_vec .* h_pattern .* kappa_o_desired);

% Solve Forward Kinematics (Force Balance)
% Unknowns: [kappa_o_1, ..., kappa_o_n, tau]
x0 = [kappa_o_desired, 10]; % Initial guess
fun = @(x) fk_solver(x, q_required, n, h_pattern, d_vec, E, I_o_vec, I_i_vec, gamma_o_vec);
sol = fsolve(fun, x0, options);

kappa_o_actual = sol(1:n);
tau_actual = sol(end);

% Convert back to centerline curvature (Eq. 9)
kappa_actual = kappa_o_actual ./ (1 + gamma_o_vec .* kappa_o_actual);

%% ================= 6. Visualization =================
% Helper to generate points
pts_desired = generate_shape(n, kappa_desired, l_segment, c, beta);
pts_actual = generate_shape(n, kappa_actual, l_segment, c, beta); % Uses fixed l_segment for comparison visual
% Note: In reality, l_j changes slightly with actuation (Eq. 10), but for visual verification 
% of curvature match, using the nominal length is usually sufficient. 
% strictly: l_j_act = h_pattern ./ (1 - gamma_o_vec .* kappa_actual);

figure('Color','w','Name','CPPR Validation');
plot(pts_desired(1,:), pts_desired(3,:), 'g--', 'LineWidth', 2, 'DisplayName', 'Desired');
hold on;
plot(pts_actual(1,:), pts_actual(3,:), 'r-o', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Model FK');
% Draw tubes (schematic)
for j=1:n
    scatter(pts_actual(1,j+1), pts_actual(3,j+1), 50, 'k.'); 
end

grid on; axis equal;
xlabel('X (mm)'); ylabel('Z (mm)');
title(sprintf('CPPR Design Validation (q = %.2f mm)', q_required));
legend('Location','best');

% Error Metrics
pos_error = vecnorm(pts_desired - pts_actual);
fprintf('Max Position Error: %.5f mm\n', max(pos_error));
fprintf('Actuation Force (Tau): %.2f N\n', tau_actual);

%% ================= Helper Functions =================

function [pts] = generate_shape(n, kappa, l_seg, c, beta)
    pts = [0;0;0];
    T = eye(4);
    for j = 1:n
        % 1. Rigid Section (c)
        T = T * [eye(3), [0;0;c]; 0 0 0 1];
        
        % 2. Curved Section (l_seg)
        kj = kappa(j);
        theta = kj * l_seg;
        if abs(kj) < 1e-8
            Tb = [eye(3), [0;0;l_seg]; 0 0 0 1];
        else
            sin_t = sin(theta); cos_t = cos(theta); inv_k = 1/kj;
            % Planar bending in X-Z plane
            Tb = [cos_t  0  sin_t  (1-cos_t)*inv_k;
                  0      1  0      0;
                 -sin_t  0  cos_t  sin_t*inv_k;
                  0      0  0      1];
        end
        T = T * Tb;
        pts(:,end+1) = T(1:3,4);
    end
end

function F = fk_solver(x, q_target, n, h, d, E, I_o, I_i, gamma_o)
    k_o = x(1:n);
    tau = x(n+1);
    
    % Eq. 7: Moment Balance at each notch
    % tau = ( E_o*I_o*k_o + E_i*I_i* (k_o / (1 - d*k_o)) ) / d
    % Rearrange to F = 0
    k_i = k_o ./ (1 - d .* k_o); % Eq. 4
    moment_imbalance = (E * I_o .* k_o + E * I_i .* k_i) ./ d - tau;
    
    % Eq. 8: Displacement Constraint
    disp_constraint = sum(d .* h .* k_o) - q_target;
    
    F = [moment_imbalance, disp_constraint];
end

function val = solve_alpha_discrepancy(g_o, target_alpha, E, OD_o, ID_o, OD_i, ID_i)
    % 1. Calculate Outer Props
    [I_o, gam_o, ~] = geometry_properties(OD_o, ID_o, g_o);
    
    % 2. Find matching Inner Props (Stiffness Match Eq. 20)
    % We need to find g_i such that I_i(g_i) == I_o
    solve_stiff = @(g) geometry_properties(OD_i, ID_i, g) - I_o;
    % Bounds check for fzero could be added
    try
        g_i = fsolve(solve_stiff, OD_i/2, optimoptions('fsolve','Display','off'));
    catch
        g_i = OD_i/2; % Fallback
    end
    
    [I_i, gam_i, ~] = geometry_properties(OD_i, ID_i, g_i);
    
    % 3. Calculate Alpha
    d = abs(gam_o) + abs(gam_i);
    EI_sum = E*(I_o + I_i);
    alpha_calc = d / EI_sum;
    
    val = alpha_calc - target_alpha;
end

function max_strain = calculate_max_strain(g, kappa, OD, ID)
    % Calculate geometric properties
    [~, gamma_mag, ~] = geometry_properties(OD, ID, g);
    
    % Strain on outer surface (Edge a) - Eq. 16
    % Distance from neutral axis to outer edge = OD/2 - gamma
    eps_a = abs(kappa) * (OD/2 - gamma_mag);
    
    % Strain on inner notch edge (Edge b) - Eq. 17
    % Distance from neutral axis to cut depth = | gamma - (OD/2 - g) |
    % Note: Neutral axis is shifted towards uncut side.
    % y_bar (gamma) is measured from center.
    % Cut is at y = OD/2 - g (from top). 
    % Distance = gamma + (OD/2 - g)
    eps_b = abs(kappa) * abs( (OD/2 - g) + gamma_mag );
    
    max_strain = max(eps_a, eps_b);
end

function [I, gamma_mag, Area] = geometry_properties(OD, ID, g)
    % Calculates Area Moment of Inertia (I) and Neutral Axis Offset (gamma)
    % for a circular tube with a chord cut of depth g.
    % Reference: Mechanics of Materials for Circular Segment / Sector
    
    R = OD/2; r = ID/2;
    
    % Check if cut goes through inner lumen
    % Usually CPPRs are deep cut, g > (OD-ID)/2
    
    % Function for properties of a solid circle segment cut at depth g
    [As, Is, ys] = circ_segment_props(R, g);
    
    % If cut is deep enough to cut inner lumen (which it always is in CPPR)
    % We treat it as: (Outer Segment) - (Inner Segment)
    % Inner cut depth g_in = g - (R - r)
    g_in = g - (R - r);
    
    if g_in <= 0
        % Not cutting inner lumen (shallow cut)
        [As_in, Is_in, ys_in] = deal(pi*r^2, pi*r^4/4, 0); % Full inner circle
    else
        [As_in, Is_in, ys_in] = circ_segment_props(r, g_in);
    end
    
    % Net Properties (Parallel Axis Theorem)
    Area = As - As_in;
    
    % Centroid from center (signed, positive away from cut)
    % moment of area / total area
    gamma_mag = (As * (-ys) - As_in * (-ys_in)) / Area;
    % Note: circ_segment_props returns centroid distance from CENTER towards cut
    % But gamma in CPPR paper is usually defined shift AWAY from cut.
    % Let's stick to magnitude here.
    gamma_mag = abs(gamma_mag);
    
    % Inertia about new neutral axis
    % I_new = I_center - A*gamma^2  <-- Wait, Parallel axis theorem:
    % I_center = I_centroid + A*d^2
    % We have I about Center.
    % I_centroid = I_center - A * y_bar^2
    
    % I calculated at center for Outer and Inner
    I_center_net = Is - Is_in;
    
    % Shift to new centroid (gamma)
    I = I_center_net - Area * gamma_mag^2;
end

function [A, I_center, y_bar] = circ_segment_props(R, g)
    % Properties of a circular area with a segment REMOVED of depth g
    % The "shape" is the remaining part.
    % But standard formulas are usually for the Segment itself.
    % Let's calculate the Full Circle and subtract the Segment (the air).
    
    if g >= 2*R; A=0; I_center=0; y_bar=0; return; end
    if g <= 0; A=pi*R^2; I_center=pi*R^4/4; y_bar=0; return; end
    
    % Properties of the REMOVED Segment (Air)
    % Height of segment = g
    % Radius = R
    % Angle alpha (half angle from center)
    d = R - g; % distance from center to chord
    alpha = acos(d/R); % radians
    
    % Area of Air Segment
    A_seg = R^2 * (alpha - sin(alpha)*cos(alpha));
    
    % Centroid of Air Segment (distance from center)
    y_seg = (2*R/3) * (sin(alpha)^3 / (alpha - sin(alpha)*cos(alpha)));
    
    % Inertia of Air Segment about Center (Ix)
    % formula: Ix = R^4/4 * (alpha - sin*cos + 2*sin^3*cos) ... complicated
    % Easier: I_x_centroid = R^4/4 * (alpha - sin alpha cos alpha) + ... 
    % Let's use simpler formulation:
    % Ix = Integral y^2 dA. 
    % Standard formula about Center:
    I_seg_center = (R^4/4) * (alpha - sin(alpha)*cos(alpha) + 2*sin(alpha)^3*cos(alpha));
    
    % Properties of Full Circle
    A_full = pi * R^2;
    I_full = pi * R^4 / 4;
    
    % Properties of Remaining Shape
    A = A_full - A_seg;
    
    % Centroid of Remaining Shape (Moment balance)
    % 0 = A * y_bar + A_seg * y_seg
    y_bar = -(A_seg * y_seg) / A; % Negative implies away from the cut
    
    % Inertia of Remaining Shape about Center
    I_center = I_full - I_seg_center;
end