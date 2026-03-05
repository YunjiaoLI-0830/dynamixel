%% CPPR Design Optimization: Rediscovering Table II
% 目的：从零开始优化 h 和 g，以拟合 Fig 11 的目标形状
clear; clc; close all;

% --- 1. 准备阶段：定义机器人基础 ---
seg_params = CPPRSegment(4.02, 3.70, 3.50, 3.20, 75e3); 
n = 10; % Table II 使用 10 个切口

% --- 2. 生成“目标数据” (The Target Curve) ---
% 我们先用"作弊"的方式，用 Table II 的参数生成目标曲线
% 假设 Table II 的参数如下 (近似值，用于生成 Target):
% h 随 index 增加 (6.02 -> 7.44)
target_h = linspace(6.02, 7.44, n); 
target_g = ones(1, n) * 2.6;         
target_c = 1.0;

% 实例化“完美机器人”并计算目标形状
robot_gt = CPPRRobot(seg_params, struct('h', target_h, 'g_o', target_g, 'c', target_c, 'g_i', zeros(1,n), 'phi_o', zeros(1,n)));
q_test = 8.0; % 假设 Fig 11 是在 q=8mm 时的形状
kappa_gt = robot_gt.solve_shape(q_test);
points_target = robot_gt.get_backbone_shape(kappa_gt); % 这就是我们要拟合的“Fig 11 曲线”

fprintf('目标曲线生成完毕。末端位置: [%.2f, %.2f, %.2f]\n', points_target(:, end));

% --- 3. 初始化优化问题 (Start from Scratch) ---
% 现在我们假装不知道 h 和 g，随便给一组初值（均匀分布）
init_h = ones(1, n) * 6.5; % 猜一个平均值
init_g = ones(1, n) * 2.0; % 猜一个浅一点的切深

% 初始机器人
robot_opt = CPPRRobot(seg_params, struct('h', init_h, 'g_o', init_g, 'c', target_c, 'g_i', zeros(1,n), 'phi_o', zeros(1,n)));

% 优化变量 x = [h_1 ... h_n, g_1 ... g_n] (共 20 个变量)
x0 = [init_h, init_g];

% 变量边界 (Bounds)
lb_h = ones(1, n) * 4.0; ub_h = ones(1, n) * 8.0; % h 限制在 4-8mm
lb_g = ones(1, n) * 0.5; ub_g = ones(1, n) * 3.0; % g 限制在 0.5-3mm
lb = [lb_h, lb_g];
ub = [ub_h, ub_g];

% --- 4. 运行优化器 (Optimization) ---
fprintf('开始参数优化，拟合目标曲线...\n');
options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxFunctionEvaluations', 2000);

% 定义代价函数: 最小化 (当前形状 - 目标形状) 的误差
cost_func = @(x) shape_fitting_cost(x, robot_opt, q_test, points_target);

[x_opt, resnorm] = lsqnonlin(cost_func, x0, lb, ub, options);

% --- 5. 解析结果与绘图 ---
h_opt = x_opt(1:n);
g_opt = x_opt(n+1:end);

fprintf('--------------------------------\n');
fprintf('优化完成!\n');
fprintf('目标 h: '); disp(target_h);
fprintf('优化 h: '); disp(h_opt);
fprintf('目标 g: '); disp(target_g);
fprintf('优化 g: '); disp(g_opt);

% 验证结果
robot_opt = robot_opt.update_design(h_opt, g_opt); % 需要你在类里把 update_design 改一下支持同时更新 h 和 g
k_final = robot_opt.solve_shape(q_test);
points_final = robot_opt.get_backbone_shape(k_final);

% 画图对比
figure('Color', 'w'); hold on; grid on; axis equal;
plot3(points_target(1,:), points_target(3,:), points_target(2,:), 'k-', 'LineWidth', 3, 'DisplayName', 'Target (Fig 11)');
plot3(points_final(1,:), points_final(3,:), points_final(2,:), 'r--', 'LineWidth', 2, 'DisplayName', 'Optimized Result');
legend;
title('Inverse Design: Fitting Fig 11 Curve');
xlabel('X'); ylabel('Z'); view(0, 90); % 侧视图

%% --- 代价函数 ---
function F = shape_fitting_cost(x, robot, q, target_pts)
    n = robot.n;
    current_h = x(1:n);
    current_g = x(n+1:end);
    
    % 更新机器人参数
    % 注意：这里需要你修改 CPPRRobot 类，添加支持更新 h 的方法
    % 或者在这里临时重新实例化 (速度较慢但代码简单)
    notches = struct('h', current_h, 'g_o', current_g, 'c', robot.c, 'g_i', robot.g_i, 'phi_o', robot.phi_o);
    temp_robot = CPPRRobot(robot.seg, notches); 
    
    try
        k = temp_robot.solve_shape(q);
        pts = temp_robot.get_backbone_shape(k);
        
        % 计算误差：所有关键点坐标差的平方
        % target_pts 是 [3 x N+1]
        diff = pts - target_pts;
        F = diff(:); % lsqnonlin 需要返回向量
    catch
        F = ones(size(target_pts(:))) * 1e5; % 惩罚
    end
end