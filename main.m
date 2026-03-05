%% CPPR Robot Simulation Main Script
clear; clc; close all;

% 1. 定义管材参数 (Table I & Caption)
% 使用 CPPRSegment 类
seg_params = CPPRSegment(...
    4.02, ...          % od_o
    4.02 - 2*0.16, ... % id_o
    3.50, ...          % od_i
    3.50 - 2*0.15, ... % id_i
    75e3 ...           % E 的值文中没有提到
);

% 2. 定义切口模式 (Table I)
notches_config.h = [2.17, 2.28, 2.24, 2.20, 2.11, 2.04, 2.16, 2.24, 2.22, 2.35];
notches_config.g_o = [3.63, 3.70, 3.68, 3.65, 3.55, 3.32, 3.62, 3.68, 3.66, 3.72];
notches_config.g_i = [3.09, 3.16, 3.14, 3.12, 3.01, 2.77, 3.09, 3.14, 3.13, 3.18];
% 前5个是 pi，后5个是 0
notches_config.phi_o = [pi, pi, pi, pi, pi, 0, 0, 0, 0, 0];
notches_config.c = 3.27;

% 3. 实例化机器人
robot = CPPRRobot(seg_params, notches_config);

% 4. 仿真与绘图
q_range = [5, 10, 15, 20]; % 驱动位移 (mm)
colors = lines(length(q_range));

figure('Color', 'w', 'Name', 'CPPR Robot Simulation (MATLAB OOP)');
hold on; grid on; axis equal;

for i = 1:length(q_range)
    q = q_range(i);
    
    % 求解形状
    kappa_sol = robot.solve_shape(q);
    
    % 获取骨架坐标
    pts = robot.get_backbone_shape(kappa_sol);
    
    % 绘图
    plot(pts(1, :), pts(3, :), 'o-', ...
        'Color', colors(i, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('q = %d mm', q));
end

xlabel('X Position (mm)');
ylabel('Z Position (mm)');
title('CPPR Robot Shape Simulation (MATLAB Refactored)');
legend('Location', 'best');
view(0, 90); % 确保 X-Z 平面视角正确