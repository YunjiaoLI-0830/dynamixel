
clear; clc; close all;
%% 1. 目标曲线生成（Target Generation）
% 目标曲线：多段常曲率圆弧拼接（piecewise constant curvature）
% 每段可定义：
%   - 段长 L_seg（mm）
%   - 弯曲角度 theta_seg（rad，正负表示方向：正=左弯，负=右弯）
%
% 例子：先左弯 π/16，再右弯 π/8，并且每段长度可以指定
% 向右弯为正，向左弯为负

segments = [
    struct('L', 20, 'theta', +15*pi/16);   % 第1段：长度60mm，向右弯 
    struct('L', 20, 'theta', +3*pi/16);   % 第2段：长度40mm，向左弯
    % struct('L', 40, 'theta', +pi/4);
%     struct('L', 60, 'theta', +pi/8);
];
% ------------------------------------------------------

% 总长度由各段长度自动得到
Total_Length = sum([segments.L]);

n_nodes = 10;  % notch 数量（最终点数为 n_nodes+1）

% 依据各段长度比例给每段分配离散段数（保证总和= n_nodes）
L_all = [segments.L];
n_each = max(1, round(n_nodes * L_all / sum(L_all)));   % 初步分配
% 修正：保证 sum(n_each) == n_nodes
while sum(n_each) > n_nodes
    [~, idx] = max(n_each);
    n_each(idx) = n_each(idx) - 1;
end
while sum(n_each) < n_nodes
    [~, idx] = max(L_all);
    n_each(idx) = n_each(idx) + 1;
end

% 使用齐次变换累积拼接每一段（在 x-z 平面内）
T = eye(4);
pts_list = [];  % 将存储 [3 x (n_nodes+1)] 的点

for k = 1:numel(segments)
    Lseg = segments(k).L;
    th   = segments(k).theta;  % 带符号角度：正左弯，负右弯

    % 该段的曲率（常曲率）：kappa = theta / L
    % 注意：theta 可以为负 => kappa 也为负 => 自动表示反方向弯曲
    kappa = th / Lseg;

    % 本段内部采样点（包含起点和终点）
    s_local = linspace(0, Lseg, n_each(k)+1);

    % 为避免相邻段重复拼接点：除第一段外，去掉每段的第一个点
    if k > 1
        s_local = s_local(2:end);
    end

    for s = s_local
        if abs(kappa) < 1e-9
            % 近似直线段
            x = 0;
            z = s;
        else
            % 常曲率圆弧的局部坐标（与机器人 forward kinematics 一致的形式）
            theta_s = kappa * s;
            x = (1 - cos(theta_s)) / kappa;
            z = sin(theta_s) / kappa;
        end

        p_local = [x; 0; z; 1];
        p_global = T * p_local;

        pts_list = [pts_list, p_global(1:3)];
    end

    % 更新当前段末端的齐次变换（用于下一段起点）
    if abs(kappa) < 1e-9
        T_end = [1,0,0,0;
                 0,1,0,0;
                 0,0,1,Lseg;
                 0,0,0,1];
    else
        theta_end = kappa * Lseg;
        s_end = sin(theta_end);
        c_end = cos(theta_end);
        inv_k = 1 / kappa;

        % 与 get_backbone_shape() 里 Tb 一致的平面弯曲变换
        T_end = [ c_end, 0,  s_end, inv_k*(1-c_end);
                  0,     1,  0,     0;
                 -s_end, 0,  c_end, inv_k*s_end;
                  0,     0,  0,     1];
    end

    T = T * T_end;
end

% 整理为目标点（3 x (n_nodes+1)）
target_points = pts_list;

% 新加入
phi_guess = estimate_phi_from_target(target_points);

fprintf('目标设定: 多段圆弧拼接，总长度 %.0f mm，总点数 %d\n', Total_Length, size(target_points,2));


%% 2. 系统初始化（System Initialization）
% A. 定义管材几何与材料参数

od_o = 4.02; id_o = 3.70;     % 外管外径 / 内径
od_i = 3.50; id_i = 3.20;     % 内管外径 / 内径
E_nitinol = 75e3;             % 镍钛合金弹性模量（MPa）

% 创建 CPPRSegment 对象，负责几何与截面属性计算
seg = CPPRSegment(od_o, id_o, od_i, id_i, E_nitinol);

% B. 初始猜测（Initial Guess）
guess_h = ones(1, n_nodes) * 2.5;   % 初始切口长度 h = 2.5 mm
guess_g = ones(1, n_nodes) * 1.5;   % 初始切口深度 g = 1.5 mm
guess_q = 1.0;                      % 初始推拉量 q = 1.0 mm

% 其他固定参数
c_fixed = 1.0;                      % 刚性段长度（mm）

% 构造 notch 参数结构体
notches_struct = struct( ...
    'h', guess_h, ...
    'g_o', guess_g, ...
    'c', c_fixed, ...
    'g_i', zeros(1, n_nodes), ...   % 内管不切割
    'phi_o', phi_guess);
%     'phi_o', zeros(1, n_nodes));    % 平面弯曲（切口方向为 0）
%       新修改

% 创建 CPPRRobot 对象
robot = CPPRRobot(seg, notches_struct);

%% 3. 优化器配置（Optimization Setup）
fprintf('启动联合优化 (h, g, q)...\n');

% 优化变量向量：
% x = [h_1 ... h_n, g_1 ... g_n, q]
x0 = [guess_h, guess_g, guess_q];

%% 变量约束（Bounds）
% 切口长度 h 的物理范围
lb_h = ones(1, n_nodes) * 2.0;
ub_h = ones(1, n_nodes) * 8.0;

% 切口深度 g 的物理范围
lb_g = ones(1, n_nodes) * 0.5;
ub_g = ones(1, n_nodes) * 3.0;

% 推拉驱动量 q 的范围
lb_q = -15.0;
% lb_q = 0;
ub_q = 15.0;

lb = [lb_h, lb_g, lb_q];
ub = [ub_h, ub_g, ub_q];

%% 优化器选项设置（lsqnonlin）
% 使用 trust-region-reflective 算法以支持边界约束
options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'trust-region-reflective', ...
    'FunctionTolerance', 1e-8, ...
    'StepTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 20000, ...
    'FiniteDifferenceType', 'central');

% 定义目标函数（残差向量）
loss_fcn = @(x) objective_function(x, robot, target_points);

% 执行优化
[x_opt, resnorm, ~, exitflag] = lsqnonlin(loss_fcn, x0, lb, ub, options);

%% 4. 结果计算与 RMS 误差分析
% 提取优化后的参数
h_final = x_opt(1:n_nodes);
g_final = x_opt(n_nodes+1 : 2*n_nodes);
q_final = x_opt(end);

% 根据优化结果重建最终形状
robot.update_design(h_final, g_final);
k_final = robot.solve_shape(q_final);
pts_final = robot.get_backbone_shape(k_final);

% --- 计算几何位置 RMS 误差 ---
% 每个点的欧氏距离误差
position_errors = sqrt(sum((pts_final - target_points).^2, 1));

% 均方根误差（RMS）
final_rms = sqrt(mean(position_errors.^2));

% 最大单点误差
max_error = max(position_errors);

fprintf('------------------------------------\n');
fprintf('optimization finished！（ExitFlag: %d）\n', exitflag);
fprintf('the acuation of q: %.4f mm\n', q_final);
fprintf('------------------------------------\n');
fprintf('position error\n');
fprintf('  RMS error:       %.5f mm\n', final_rms);
fprintf('  max error:    %.5f mm\n', max_error);
fprintf('------------------------------------\n');

%% 5. 结果可视化
figure('Color', 'w', 'Position', [100, 100, 1000, 400]);

% 子图 1：目标曲线与优化结果对比
subplot(1, 2, 1); hold on; grid on; axis equal;
plot3(target_points(1,:), target_points(3,:), target_points(2,:), ...
      'k-', 'LineWidth', 2, 'DisplayName', 'Target');
plot3(pts_final(1,:), pts_final(3,:), pts_final(2,:), ...
      'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', ...
      'DisplayName', 'Optimized');
title(sprintf(' (RMS: %.4f mm)', final_rms));
xlabel('X'); zlabel('Z'); view(0, 90);
legend('Location','best');

% 子图 2：优化后的设计参数分布
subplot(1, 2, 2); hold on; grid on; box on;
yyaxis left
plot(1:n_nodes, h_final, 's-', 'LineWidth', 2);
ylabel('notches length h [mm]');
yyaxis right
plot(1:n_nodes, g_final, 'd-', 'LineWidth', 2);
ylabel('notches depth g [mm]');
xlabel('notches number');
title('optimization parameters');

%% 6. 目标函数定义
function F = objective_function(x, robot, target_pts)
    % 从优化变量中解析 h、g、q
    n = robot.n;
    h_c = x(1:n);
    g_c = x(n+1 : 2*n);
    q_c = x(end);
    
    % 更新机器人几何设计参数
    robot.update_design(h_c, g_c);
    
    try
        % 正运动学求解
        k = robot.solve_shape(q_c);
        pts = robot.get_backbone_shape(k);
        
        % 几何位置残差（用于形状拟合）
        diff_geom = (pts - target_pts);
        w_smooth = 0.05;
        diff_h = diff(h_c) * w_smooth;
        diff_g = diff(g_c) * w_smooth;
        
%         拼接所有残差项
        F = [diff_geom(:); diff_h(:); diff_g(:)];
        
    catch
        % 若正运动学求解失败，返回大残差以惩罚该解
        F = ones(numel(target_pts) + 2*(n-1), 1) * 1e5;
    end
end

function phi_o = estimate_phi_from_target(target_points)
    % 根据离散点的转向（x-z 平面）估计每段弯曲方向
    % 输出：1 x n，其中每个元素为 0（左弯）或 pi（右弯）
    %
    % 做法：计算相邻段的转角符号（2D 叉积符号）

    x = target_points(1,:);
    z = target_points(3,:);
    N = size(target_points,2);
    n = N - 1;

    % 切向量
    dx = diff(x);
    dz = diff(z);

    % 默认全 0
    phi_o = zeros(1,n);

    % 对每个“节点段”估计转向：用 v_i 和 v_{i+1} 的叉积符号
    for i = 1:n
        if i == 1
            v1 = [dx(i); dz(i)];
            v2 = [dx(i); dz(i)];
        elseif i == n
            v1 = [dx(i-1); dz(i-1)];
            v2 = [dx(i);   dz(i)];
        else
            v1 = [dx(i-1); dz(i-1)];
            v2 = [dx(i);   dz(i)];
        end

        cross2 = v1(1)*v2(2) - v1(2)*v2(1); % 2D 叉积（标量）
        if cross2 < 0
            phi_o(i) = pi;   % 右弯
        else
            phi_o(i) = 0;    % 左弯或直线
        end
    end

    % 平滑一下方向翻转：避免单点噪声导致 0/pi 抖动（可选但很有用）
    for i = 2:n-1
        if phi_o(i-1)==phi_o(i+1) && phi_o(i)~=phi_o(i-1)
            phi_o(i) = phi_o(i-1);
        end
    end
end


