%% CPPR validation
clear; clc; close all;

%% ================= Inverse Design  =================

n = 10;                 
E = 75e3;               
OD_o = 4.02;            % 目标外管直径 (mm)
OD_i = 3.50;            % 目标外管直径 (mm)
epsilon_max = 0.02;     % 材料极限 2% for Nitinol
FS = 1.5;               % 安全系数
epsilon_FS = epsilon_max / FS;  % 设计应变极限 
l_segment = 5.0;         % 设定每节的目标弧长

% sampling of the S shape
kappa_desired = [0.01, 0.02, 0.03, 0.02, 0.01, -0.01, -0.02, -0.03, -0.02, -0.01];

%% 2. find the Benchmark Notch
%  maximum of kappa
[max_kappa, m_idx] = max(abs(kappa_desired));

% maximum of kappa
% Eq. 22

g_initial = OD_o * 0.8;
options = optimoptions('fsolve', 'Display', 'off');
% strain limit constraint
% 需要改动
% ！！！
% 计算偏移量(g-OD_-/2)*0.4不准确
g_bench = fsolve(@(g) (abs(max_kappa)*(OD_o/2 - (g-OD_o/2)*0.4) - epsilon_FS), g_initial, options);

%% 3. calculate the depth scaling
g_pattern = zeros(1, n);
phi_pattern = zeros(1, n);

for j = 1:n
    % Eq.19
    g_pattern(j) = g_bench * abs(kappa_desired(j) / max_kappa);
    
    % the oriention of notches
    % !!
    if kappa_desired(j) >= 0
        phi_pattern(j) = pi;
    else
        phi_pattern(j) = 0;
    end
end

%% 4. calculate h
h_pattern = zeros(1, n);
gamma_o_vec = zeros(1, n);

for j = 1:n
    % Eq.10
    % 带方向！
    gamma_o_val = (g_pattern(j) - OD_o/2) * 0.4 * cos(phi_pattern(j));
    gamma_o_vec(j) = gamma_o_val;
    
    % notch length Eq.10
    % 与原文不同！
    h_pattern(j) = l_segment * (1 - gamma_o_val * kappa_desired(j));
end

%% ================= Forward Verification =================
% 简化Eq.20
g_i_pattern = g_pattern * 0.85; 

% calculate offset of netral line and I
% cite[24,27]
gamma_o_fk = (g_pattern - OD_o/2) * 0.4 .* cos(phi_pattern); 
gamma_i_fk = (g_i_pattern - OD_i/2) * 0.4 .* cos(phi_pattern - pi); 

% Eq.3: 力臂 d
d_vec = gamma_i_fk - gamma_o_fk; 

I_o_vec = 0.5 * pi * (OD_o/2)^4 .* (1 - (g_pattern/OD_o)); 
I_i_vec = 0.5 * pi * (OD_i/2)^4 .* (1 - (g_i_pattern/OD_i));

% Eq.8
% 
q_required = sum(d_vec .* h_pattern .* kappa_desired);

% 求解力平衡方程
x0 = [kappa_desired, 1]; % 初始猜测
fun = @(x) solver_helper(x, q_required, n, h_pattern, d_vec, E, I_o_vec, I_i_vec);
sol = fsolve(fun, x0, options);
kappa_actual = sol(1:n); % 反解出来的实际曲率

%% ================= Visualization  =================

% Desired Curve
pts_desired = [0;0;0]; T = eye(4);
for j = 1:n
    kj = kappa_desired(j);
    lj = l_segment; % 目标弧长是固定的 5.0mm
    theta = kj*lj;
    if abs(kj)<1e-8, Tb=eye(4); Tb(3,4)=lj; else
        Tb = [cos(theta) 0 sin(theta) (1-cos(theta))/kj; 0 1 0 0; -sin(theta) 0 cos(theta) sin(theta)/kj; 0 0 0 1];
    end
    T = T * Tb * [1 0 0 0; 0 1 0 0; 0 0 1 3.27; 0 0 0 1]; % 加上刚性段 c=3.27
    pts_desired(:,end+1) = T(1:3,4);
end

% Solve Curve
pts_actual = [0;0;0]; T = eye(4);
for j = 1:n
    kj_o = kappa_actual(j);
    gamma_j = gamma_o_fk(j);
    
    % Eq.9 & Eq.10 
    kj = kj_o / (1 + gamma_j * kj_o);
    lj = h_pattern(j) / (1 - gamma_j * kj); % 长度修正
    
    theta = kj*lj;
    if abs(kj)<1e-8, Tb=eye(4); Tb(3,4)=lj; else
        Tb = [cos(theta) 0 sin(theta) (1-cos(theta))/kj; 0 1 0 0; -sin(theta) 0 cos(theta) sin(theta)/kj; 0 0 0 1];
    end
    T = T * Tb * [1 0 0 0; 0 1 0 0; 0 0 1 3.27; 0 0 0 1];
    pts_actual(:,end+1) = T(1:3,4);
end

% plot
figure('Color','w','Name','Closed-Loop Validation');
plot(pts_desired(1,:), pts_desired(3,:), 'g--', 'LineWidth', 3, 'DisplayName', 'Desired Curve (Input)');
hold on;
plot(pts_actual(1,:), pts_actual(3,:), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Actual Solved Curve (FK)');
grid on; axis equal;
xlabel('X Position (mm)'); ylabel('Z Position (mm)');
title('Inverse Design -> Forward Kinematics Verification');
legend('Location','best');

% print
fprintf('validation result:\n');
fprintf('q_required: %.2f mm\n', q_required);
fprintf('the maxmum of error: %.4f mm\n', max(vecnorm(pts_desired - pts_actual)));

%% ================= 辅助函数 =================
function F = solver_helper(x, q_target, n, h, d, E, I_o, I_i)
    kappa_o = x(1:n);
    tau = x(n+1);
    F = zeros(1, n+1);
    for j = 1:n
        % Eq.7 静力学平衡
        F(j) = (E*I_o(j)*kappa_o(j) + E*I_i(j)*(kappa_o(j)/(1 - d(j)*kappa_o(j)))) / d(j) - tau;
    end
    % Eq.8 位移约束
    F(n+1) = sum(d .* h .* kappa_o) - q_target;
end