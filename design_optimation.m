%% CPPR 逆向设计优化演示 (Inverse Design Optimization)
clear; clc; close all;

% 1. 初始设置 (基础参数)
% 定义管材常数 (Nitinol)
seg_params = CPPRSegment(4.02, 3.70, 3.50, 3.20, 75e3); 

% 初始猜测的设计参数 (均匀切深)
n_notches = 5;
initial_g = ones(1, n_notches) * 0.5 * (4.02/2); % 初始切深为半径的一半
notches_init.h = ones(1, n_notches) * 2;   % 切口长 2mm
notches_init.c = 1;                        % 间距 1mm
notches_init.g_o = initial_g; 
notches_init.g_i = zeros(1, n_notches);    % 假设内管不切或固定
notches_init.phi_o = zeros(1, n_notches);  % 全部朝向 0度

% 实例化机器人对象
robot = CPPRRobot(seg_params, notches_init);

% 2. 定义优化目标 (Target)
% 假设我们要让机器人到达一个特定的空间点
target_point = [15; 0; 85]; % (x, y, z) mm

% 3. 定义优化变量 x
% 变量 x 包括两部分： [ 设计参数(g_o)  |  控制参数(q) ]
% x = [g1, g2, g3, g4, g5, q] (共 n+1 个变量)
x0 = [initial_g, 0]; % 初值：当前切深 + 0驱动

% 4. 设置边界 (Bounds)
% 切深 g 的范围: 0.1mm 到 3.0mm (不能切断管子)
lb_g = ones(1, n_notches) * 0.1;
ub_g = ones(1, n_notches) * 3.5;

% 驱动 q 的范围: -10mm 到 10mm
lb_q = -10;
ub_q = 10;

lb = [lb_g, lb_q];
ub = [ub_g, ub_q];

% 5. 运行优化 (fmincon)
fprintf('开始设计优化...\n');
options = optimoptions('fmincon', ...
    'Display', 'iter', ...       % 显示过程
    'Algorithm', 'sqp', ...      % SQP 算法适合这种中小规模非线性问题
    'MaxFunctionEvaluations', 1000);

% 定义目标函数句柄
cost_func = @(x) objective_function(x, robot, target_point);

[x_opt, fval] = fmincon(cost_func, x0, [],[],[],[], lb, ub, [], options);

% 6. 解析结果与可视化
g_opt = x_opt(1:n_notches);
q_opt = x_opt(end);

fprintf('--------------------------------\n');
fprintf('优化完成!\n');
fprintf('目标点: [%.1f, %.1f, %.1f]\n', target_point);
fprintf('优化后驱动 q: %.4f mm\n', q_opt);
fprintf('优化后切深 g: \n'); disp(g_opt);

% --- 画图对比 ---
% 1. 初始设计 (Initial Design)
robot.update_design(initial_g);
k_init = robot.solve_shape(0); % 初始状态(无驱动)
pts_init = robot.get_backbone_shape(k_init);

% 2. 优化后设计 (Optimized Design)
robot.update_design(g_opt);      % 更新为最优切深
k_opt = robot.solve_shape(q_opt);% 使用最优驱动
pts_opt = robot.get_backbone_shape(k_opt);

figure('Color', 'w');
plot3(pts_init(1,:), pts_init(2,:), pts_init(3,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Initial (Straight)'); hold on;
plot3(pts_opt(1,:), pts_opt(2,:), pts_opt(3,:), 'r-', 'LineWidth', 3, 'DisplayName', 'Optimized Design');
plot3(target_point(1), target_point(2), target_point(3), 'b*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Target');
grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
legend;
title('Inverse Design Optimization Result');
view(0, 0); % 侧视图

%% --- 目标函数 (Cost Function) ---
function cost = objective_function(x, robot, target_point)
    % 1. 提取变量
    n = robot.n;
    current_g = x(1:n);
    current_q = x(n+1);
    
    % 2. 更新机器人几何 (Design Update)
    robot = robot.update_design(current_g);
    
    % 3. 运行正运动学 (Physics Solver)
    try
        kappa = robot.solve_shape(current_q);
    catch
        cost = 1e6; % 如果求解失败，给一个巨大的惩罚值
        return;
    end
    
    % 4. 获取末端位置
    points = robot.get_backbone_shape(kappa);
    tip_pos = points(:, end);
    
    % 5. 计算代价
    % 项1: 距离误差 (主要目标)
    dist_error = norm(tip_pos - target_point);
    
    % 项2: 正则化 (Regularization) - 让切深变化平滑
    % 避免 g = [0.1, 3.0, 0.1, 3.0] 这种锯齿状设计，不仅难加工而且应力集中
    smoothness = sum(diff(current_g).^2); 
    
    % 总代价 = 距离^2 + 权重 * 平滑度
    cost = dist_error^2 + 10 * smoothness;
end