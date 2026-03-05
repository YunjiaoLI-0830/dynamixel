%% 4.6 (按论文第二章 II-E) 用 sequential PCC 只拟合 kappa_ideal，并与模型 k_final 对比
% 需要：target_points (3x(N)), h_final (1xn), c_fixed, 以及 k_final (1xn)
% 说明：这里把“理想曲率”定义为：在给定(h_final, c_fixed)的离散结构下，
%      用 II-E 的逐段单标量优化，使每一步的下一点尽量贴合 target_points。

% --- 用 II-E sequential 方法估计 ideal kappa ---
[kappa_ideal_seq, pts_seq_fit] = estimate_kappa_sequential_IIE(target_points, h_final, c_fixed);

% --- 确保 model curvature 和 ideal curvature 维度一致 ---
kappa_model = k_final(:).';                 % 1 x n
kappa_ideal = kappa_ideal_seq(:).';         % 1 x n

% --- 曲率误差 ---
dk = kappa_model - kappa_ideal;
rms_k = sqrt(mean(dk.^2));
max_k = max(abs(dk));

fprintf('------------------------------------\n');
fprintf('II-E sequential PCC curvature error\n');
fprintf('  RMS(dkappa): %.6f [1/mm]\n', rms_k);
fprintf('  MAX(|dkappa|): %.6f [1/mm]\n', max_k);
fprintf('------------------------------------\n');

% --- 画图：曲率对比 + 误差 + (可选) sequential 拟合轨迹对比 ---
figure('Color','w','Position',[120,120,1200,420]);

subplot(1,3,1); hold on; grid on; box on;
plot(1:numel(kappa_ideal), kappa_ideal, 'k-', 'LineWidth', 2);
plot(1:numel(kappa_model), kappa_model, 'r--', 'LineWidth', 2);
xlabel('notch index j');
ylabel('\kappa_j [1/mm]');
title('Curvature: Ideal (II-E) vs Model');
legend('Ideal (sequential II-E)','Model (robot.solve\_shape)','Location','best');

subplot(1,3,2); hold on; grid on; box on;
plot(1:numel(dk), dk, 'b-', 'LineWidth', 2);
yline(0,'k-');
xlabel('notch index j');
ylabel('\Delta\kappa_j [1/mm]');
title(sprintf('Curvature error (RMS=%.6f)', rms_k));

subplot(1,3,3); hold on; grid on; axis equal; box on;
% x-z 平面画轨迹（你的 plot3 用的是 (x,z,y) 的排列；这里直接画 x-z 更直观）
plot(target_points(1,:), target_points(3,:), 'k-', 'LineWidth', 2);
plot(pts_seq_fit(1,:),  pts_seq_fit(3,:),  'g--', 'LineWidth', 2);
% plot(pts_final(1,:),    pts_final(3,:),    'r-', 'LineWidth', 1.5);
xlabel('x [mm]'); ylabel('z [mm]');
title('Backbone: target vs sequential-fit vs model');
legend('Target','Sequential fit (II-E)','Model backbone','Location','best');

%% -------------------- 函数：按 II-E sequential 解每段一个 kappa --------------------
function [kappa_seq, pts_fit] = estimate_kappa_sequential_IIE(target_pts, h_vec, c)
    % target_pts: 3 x (n+1) 离散目标点
    % h_vec: 1 x n notch 弯曲段长度
    % c: 刚性段长度（每段后接 rigid segment）
    %
    % 输出：
    % kappa_seq: 1 x n，每段常曲率（II-E sequential 拟合得到）
    % pts_fit:   3 x (n+1)，用拟合的 kappa 串联得到的点列（用于可视化）

    n = numel(h_vec);
    assert(size(target_pts,2) == n+1, 'target_pts 点数必须为 n+1');

    T = eye(4);
    pts_fit = zeros(3, n+1);
    p0 = T * [0;0;0;1];
    pts_fit(:,1) = p0(1:3);

    kappa_seq = zeros(1,n);

    % 给一个合理的曲率搜索上界：假设每段最大转角不超过 pi
    h_min = max(1e-6, min(h_vec));
    kmax_global = pi / h_min;   % [1/mm]

    for j = 1:n
        p_target_next = target_pts(:, j+1);

        hj = h_vec(j);

        % 用标量优化：找 kappa_j 使预测的下一点最接近目标下一点
        % 这里用 fminbnd 更稳，不要求 fzero 的符号变号
        kmax = kmax_global;
        cost = @(kappa) step_point_cost(T, kappa, hj, c, p_target_next);

        % 如果你希望更窄/更物理的曲率范围，可把 kmax 改小
        kappa_j = fminbnd(cost, -kmax, +kmax);

        kappa_seq(j) = kappa_j;

        % 更新位姿：T <- T * T_bend(kappa,h) * T_rigid(c)
        T = T * T_bend_planar(kappa_j, hj) * T_rigid_z(c);

        p_next = T * [0;0;0;1];
        pts_fit(:, j+1) = p_next(1:3);
    end
end

function J = step_point_cost(Tprev, kappa, h, c, p_target_next)
    % 用“弯曲段 + 刚段”预测下一点，然后与目标下一点做欧氏距离平方
    Tnext = Tprev * T_bend_planar(kappa, h) * T_rigid_z(c);
    p_pred = Tnext * [0;0;0;1];
    d = p_pred(1:3) - p_target_next(:);
    J = d.'*d;
end

function Tb = T_bend_planar(kappa, s)
    % 平面常曲率弯曲（x-z 平面）齐次变换
    % 与你脚本里 T_end 的形式一致：绕 y 轴转，伴随圆弧位移
    if abs(kappa) < 1e-12
        Tb = [1,0,0,0;
              0,1,0,0;
              0,0,1,s;
              0,0,0,1];
        return;
    end

    th = kappa * s;
    c = cos(th);  si = sin(th);
    invk = 1/kappa;

    Tb = [ c,  0,  si, invk*(1-c);
           0,  1,  0,  0;
          -si, 0,  c,  invk*si;
           0,  0,  0,  1];
end

function Tr = T_rigid_z(c)
    % 刚性直段：沿局部 z 方向平移 c
    Tr = [1,0,0,0;
          0,1,0,0;
          0,0,1,c;
          0,0,0,1];
end