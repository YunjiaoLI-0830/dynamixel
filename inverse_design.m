%% CPPR Inverse Design tro Oliver 2022
clear; clc;

%% 1. input：desired shape & meterial
n = 10;                 % 计划切口数量
OD = 4.02;              % 目标外管直径 (mm)
ID = 3.50;              % 目标外管内径 (mm)
epsilon_max = 0.02;     % 材料应变极限 (2% for Nitinol) 
FS = 1.5;               % 安全系数
epsilon_FS = epsilon_max / FS; % 设计应变极限 

% S shape
% sampling of the S shape
kappa_desired = [0.01, 0.02, 0.03, 0.02, 0.01, -0.01, -0.02, -0.03, -0.02, -0.01]; 

%% 2. find the Benchmark Notch
%  maximum of kappa
[max_kappa, m_idx] = max(abs(kappa_desired)); 

% maximum of kappa
% 应变公式;Eq. 22

g_initial = OD * 0.8; 
options = optimoptions('fsolve', 'Display', 'none');
g_benchmark = fsolve(@(g) strain_error(g, max_kappa, OD, ID, epsilon_FS), g_initial, options);

%% 3. calculate the Depth Scaling

g_pattern = zeros(1, n);
for j = 1:n
    % Eq.19
    target_stiffness_ratio = abs(kappa_desired(j) / max_kappa);
    g_pattern(j) = g_benchmark * target_stiffness_ratio; 
end

%% 4. calculate h
l_segment = 5.0; 
h_pattern = zeros(1, n);
for j = 1:n
    % notch length Eq.10
    gamma_val = calculate_gamma(g_pattern(j), OD, ID);
    h_pattern(j) = l_segment * (1 - gamma_val * kappa_desired(j)); % [cite: 267]
end

%% 
fprintf('inverse design completed！\n');
fprintf('suggusted benchmark of notches depth is: %.2f mm\n', g_benchmark);
fprintf('parameters of g_pattern is %.2f mm\n', g_pattern);
fprintf('parameters of h_pattern is %.2f mm\n', h_pattern);

%% strain limit constraint
% Eq.17
function err = strain_error(g, kappa, OD, ID, target_epsilon)
    gamma = calculate_gamma(g, OD, ID);
    % Eq.16
    current_epsilon = abs(kappa) * (OD/2 - gamma); 
    err = current_epsilon - target_epsilon;
end

function gamma = calculate_gamma(g, OD, ID)
    % modify! cite[24,27]
    gamma = (g - OD/2) * 0.45; 
end