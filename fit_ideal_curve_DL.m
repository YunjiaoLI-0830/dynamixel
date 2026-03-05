%% fit_ideal_curve_DL.m
% 用深度学习(监督学习)替代 lsqnonlin：
% 1) 随机采样 (h,g,q) -> 用你的 CPPRRobot 正运动学生成 shape points
% 2) 训练回归网络： target_points -> (h,g,q)
% 3) 推理时：网络给出 (h,g,q) 作为“直接解”或作为后续数值优化 warm-start
%
% 依赖：Deep Learning Toolbox
% 可选：Statistics and Machine Learning Toolbox (用于 mapminmax 也可手写)

clear; clc; close all;

%% 0) 目标曲线生成（复用你原来的代码段）
segments = [
    struct('L', 60, 'theta', -pi/3);
    struct('L', 60, 'theta', +pi/2);
    struct('L', 60, 'theta', +pi/5);
];

Total_Length = sum([segments.L]);
n_nodes = 20;  % notch 数量（最终点数为 n_nodes+1）

L_all = [segments.L];
n_each = max(1, round(n_nodes * L_all / sum(L_all)));
while sum(n_each) > n_nodes
    [~, idx] = max(n_each);
    n_each(idx) = n_each(idx) - 1;
end
while sum(n_each) < n_nodes
    [~, idx] = max(L_all);
    n_each(idx) = n_each(idx) + 1;
end

T = eye(4);
pts_list = [];
for k = 1:numel(segments)
    Lseg = segments(k).L;
    th   = segments(k).theta;
    kappa = th / Lseg;
    s_local = linspace(0, Lseg, n_each(k)+1);
    if k > 1
        s_local = s_local(2:end);
    end
    for s = s_local
        if abs(kappa) < 1e-9
            x = 0; z = s;
        else
            theta_s = kappa * s;
            x = (1 - cos(theta_s)) / kappa;
            z = sin(theta_s) / kappa;
        end
        p_global = T * [x;0;z;1];
        pts_list = [pts_list, p_global(1:3)];
    end
    if abs(kappa) < 1e-9
        T_end = [1,0,0,0; 0,1,0,0; 0,0,1,Lseg; 0,0,0,1];
    else
        theta_end = kappa * Lseg;
        s_end = sin(theta_end); c_end = cos(theta_end);
        inv_k = 1 / kappa;
        T_end = [ c_end, 0,  s_end, inv_k*(1-c_end);
                  0,     1,  0,     0;
                 -s_end, 0,  c_end, inv_k*s_end;
                  0,     0,  0,     1];
    end
    T = T * T_end;
end
target_points = pts_list; % 3 x (n_nodes+1)
phi_guess = estimate_phi_from_target(target_points);

fprintf('目标设定: 总长度 %.0f mm, 点数 %d\n', Total_Length, size(target_points,2));

%% 1) 系统初始化（与你原来一致）
od_o = 4.02; id_o = 3.70;
od_i = 3.50; id_i = 3.20;
E_nitinol = 75e3; % MPa
seg = CPPRSegment(od_o, id_o, od_i, id_i, E_nitinol);

c_fixed = 1.0;

% 变量约束（与原 lsqnonlin 一致）
lb_h = ones(1, n_nodes) * 2.0;
ub_h = ones(1, n_nodes) * 8.0;
lb_g = ones(1, n_nodes) * 0.5;
ub_g = ones(1, n_nodes) * 3.0;
lb_q = -15.0;
ub_q = 15.0;

%% 2) 数据集生成
% 你可以先用小一点的 N_samples 试运行，确认流程没问题，再加大。
N_samples = 6000;           % 数据量：越大越准，但生成越慢
max_tries_per_sample = 20;  % fsolve 失败时重试次数

notches_struct = struct('h', ones(1,n_nodes)*2.5, ...
    'g_o', ones(1,n_nodes)*1.5, ...
    'c', c_fixed, ...
    'g_i', zeros(1,n_nodes), ...
    'phi_o', phi_guess);
robot = CPPRRobot(seg, notches_struct);

fprintf('开始生成训练数据 N=%d ...\n', N_samples);
[X, Y] = generate_dataset(robot, n_nodes, target_points, ...
    lb_h, ub_h, lb_g, ub_g, lb_q, ub_q, N_samples, max_tries_per_sample);
% X: (3*(n_nodes+1)) x N  (把 target_points 展平)
% Y: (2*n_nodes+1) x N    (h,g,q)

%% 3) 归一化
% 对输入点云和输出参数做标准化（z-score）更稳定
[Xn, Xmu, Xsig] = zscore_columns(X);
[Yn, Ymu, Ysig] = zscore_columns(Y);

%% 4) 定义网络（回归）
inputSize = size(Xn,1);
outputSize = size(Yn,1);

layers = [
    featureInputLayer(inputSize, 'Normalization','none', 'Name','in')
    fullyConnectedLayer(512, 'Name','fc1')
    reluLayer('Name','relu1')
    fullyConnectedLayer(512, 'Name','fc2')
    reluLayer('Name','relu2')
    fullyConnectedLayer(256, 'Name','fc3')
    reluLayer('Name','relu3')
    fullyConnectedLayer(outputSize, 'Name','out')
    regressionLayer('Name','reg')
];

% 划分训练/验证
rng(0);
N = size(Xn,2);
idx = randperm(N);
nTrain = round(0.9*N);
tr = idx(1:nTrain);
va = idx(nTrain+1:end);

XTrain = Xn(:,tr)';  YTrain = Yn(:,tr)';
XVal   = Xn(:,va)';  YVal   = Yn(:,va)';

opts = trainingOptions('adam', ...
    'MaxEpochs', 80, ...
    'MiniBatchSize', 256, ...
    'InitialLearnRate', 1e-3, ...
    'Shuffle','every-epoch', ...
    'ValidationData', {XVal, YVal}, ...
    'ValidationFrequency', 50, ...
    'Verbose', true, ...
    'Plots','training-progress');

net = trainNetwork(XTrain, YTrain, layers, opts);

%% 5) 用训练好的网络直接预测 (h,g,q)
Xtest = reshape(target_points, [], 1);           % 3*(n_nodes+1) x 1
Xtestn = (Xtest - Xmu) ./ Xsig;
Ypredn = predict(net, Xtestn');                  % 1 x (2n+1)
Ypred  = (Ypredn' .* Ysig) + Ymu;                % 还原尺度

h_pred = Ypred(1:n_nodes).';
g_pred = Ypred(n_nodes+1:2*n_nodes).';
q_pred = Ypred(end);

% 截断到物理范围（网络输出可能略越界）
h_pred = min(max(h_pred, lb_h), ub_h);
g_pred = min(max(g_pred, lb_g), ub_g);
q_pred = min(max(q_pred, lb_q), ub_q);

%% 6) 评价预测效果：正运动学 -> RMS
robot.update_design(h_pred, g_pred);
k_pred = robot.solve_shape(q_pred);
pts_pred = robot.get_backbone_shape(k_pred);

pos_err = sqrt(sum((pts_pred - target_points).^2, 1));
rms_pred = sqrt(mean(pos_err.^2));
fprintf('\nDL 直接预测 RMS = %.5f mm, q=%.3f\n', rms_pred, q_pred);

%% 7) （可选）DL 预测作为 warm-start，再跑一次 lsqnonlin 精修
do_refine = true;
if do_refine
    x0 = [h_pred, g_pred, q_pred];
    lb = [lb_h, lb_g, lb_q];
    ub = [ub_h, ub_g, ub_q];
    options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'Algorithm', 'trust-region-reflective', ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxFunctionEvaluations', 8000, ...
        'FiniteDifferenceType', 'central');
    loss_fcn = @(x) objective_function(x, robot, target_points);
    [x_opt, ~] = lsqnonlin(loss_fcn, x0, lb, ub, options);
    h_final = x_opt(1:n_nodes);
    g_final = x_opt(n_nodes+1:2*n_nodes);
    q_final = x_opt(end);
    robot.update_design(h_final, g_final);
    k_final = robot.solve_shape(q_final);
    pts_final = robot.get_backbone_shape(k_final);
    pos_err2 = sqrt(sum((pts_final - target_points).^2, 1));
    rms2 = sqrt(mean(pos_err2.^2));
    fprintf('DL warm-start + lsqnonlin 精修 RMS = %.5f mm (q=%.3f)\n', rms2, q_final);
else
    h_final = h_pred; g_final = g_pred; q_final = q_pred;
    pts_final = pts_pred;
    rms2 = rms_pred;
end

%% 8) 可视化
figure('Color','w','Position',[100,100,1000,400]);
subplot(1,2,1); hold on; grid on; axis equal;
plot3(target_points(1,:), target_points(3,:), target_points(2,:), ...
    'k-', 'LineWidth', 2, 'DisplayName','Target');
plot3(pts_final(1,:), pts_final(3,:), pts_final(2,:), ...
    'r-o', 'LineWidth', 1.5, 'MarkerFaceColor','r', 'DisplayName','Pred/Refined');
title(sprintf('DL/RL surrogate (RMS %.4f mm)', rms2));
xlabel('X'); zlabel('Z'); view(0,90);
legend('Location','best');

subplot(1,2,2); hold on; grid on; box on;
yyaxis left
plot(1:n_nodes, h_final, 's-', 'LineWidth', 2);
ylabel('notches length h [mm]');
yyaxis right
plot(1:n_nodes, g_final, 'd-', 'LineWidth', 2);
ylabel('notches depth g [mm]');
xlabel('notches number');
title('predicted/refined parameters');

%% ----------------- 下面是局部函数 -----------------
function [X, Y] = generate_dataset(robot, n_nodes, target_points, ...
    lb_h, ub_h, lb_g, ub_g, lb_q, ub_q, N_samples, max_tries)

    inDim  = numel(target_points);
    outDim = 2*n_nodes + 1;
    X = zeros(inDim,  N_samples, 'single');
    Y = zeros(outDim, N_samples, 'single');

    % 为了更“覆盖”目标附近的空间，可以把随机采样做成：
    % - 一半完全随机
    % - 一半在某个基准附近加噪声（可选）
    for i = 1:N_samples
        ok = false;
        for t = 1:max_tries
            h = lb_h + (ub_h - lb_h) .* rand(1,n_nodes);
            g = lb_g + (ub_g - lb_g) .* rand(1,n_nodes);
            q = lb_q + (ub_q - lb_q) .* rand(1,1);

            robot.update_design(h, g);
            try
                k = robot.solve_shape(q);
                pts = robot.get_backbone_shape(k);

                % 输入：把“shape points”展开成向量
                X(:,i) = single(reshape(pts, [], 1));
                % 输出：对应的 (h,g,q)
                Y(:,i) = single([h, g, q])';
                ok = true;
                break;
            catch
                ok = false;
            end
        end
        if ~ok
            % 极少数情况下持续失败：用 target_points + 随机输出填充，避免中断
            X(:,i) = single(reshape(target_points, [], 1));
            Y(:,i) = single([lb_h, lb_g, 0])';
        end
        if mod(i, 250) == 0
            fprintf('  数据生成进度: %d / %d\n', i, N_samples);
        end
    end
end

function [Xn, mu, sig] = zscore_columns(X)
    mu = mean(X, 2);
    sig = std(X, 0, 2);
    sig(sig < 1e-12) = 1;
    Xn = (X - mu) ./ sig;
end

function F = objective_function(x, robot, target_pts)
    n = robot.n;
    h_c = x(1:n);
    g_c = x(n+1:2*n);
    q_c = x(end);
    robot.update_design(h_c, g_c);
    try
        k = robot.solve_shape(q_c);
        pts = robot.get_backbone_shape(k);
        diff_geom = (pts - target_pts);
        w_smooth = 0.05;
        diff_h = diff(h_c) * w_smooth;
        diff_g = diff(g_c) * w_smooth;
        F = [diff_geom(:); diff_h(:); diff_g(:)];
    catch
        F = ones(numel(target_pts) + 2*(n-1), 1) * 1e5;
    end
end

function phi_o = estimate_phi_from_target(target_points)
    x = target_points(1,:);
    z = target_points(3,:);
    N = size(target_points,2);
    n = N - 1;
    dx = diff(x);
    dz = diff(z);
    phi_o = zeros(1,n);
    for i = 1:n
        if i == 1
            v1 = [dx(i); dz(i)]; v2 = [dx(i); dz(i)];
        elseif i == n
            v1 = [dx(i-1); dz(i-1)]; v2 = [dx(i); dz(i)];
        else
            v1 = [dx(i-1); dz(i-1)]; v2 = [dx(i); dz(i)];
        end
        cross2 = v1(1)*v2(2) - v1(2)*v2(1);
        if cross2 < 0
            phi_o(i) = pi;
        else
            phi_o(i) = 0;
        end
    end
    for i = 2:n-1
        if phi_o(i-1)==phi_o(i+1) && phi_o(i)~=phi_o(i-1)
            phi_o(i) = phi_o(i-1);
        end
    end
end