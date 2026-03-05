%% Forward Kinematics tro 2022 Oliver
clear; clc; close all;

%% Table 1
n = 10;                 % number of notches
E = 75e3;               % 镍钛合金弹性模量
OD_o = 4.02;            % 外管外径 (mm) 
OD_i = 3.50;            % 内管外径 (mm) 
c = 3.27;               % 切口间距 (mm) 

% 每个切口的具体几何数据
h = [2.17, 2.28, 2.24, 2.20, 2.11, 2.04, 2.16, 2.24, 2.22, 2.35]; % 切口长度 (mm)
g_o = [3.63, 3.70, 3.68, 3.65, 3.55, 3.32, 3.62, 3.68, 3.66, 3.72]; % 外管切深 (mm)
g_i = [3.09, 3.16, 3.14, 3.12, 3.01, 2.77, 3.09, 3.14, 3.13, 3.18]; % 内管切深 (mm)
phi_o = [pi, pi, pi, pi, pi, 0, 0, 0, 0, 0]; % 外管切口方向 

% calculate offset of netral line and I
% cite[24,27]
gamma_o = (g_o - OD_o/2) * 0.4 .* cos(phi_o); 
gamma_i = (g_i - OD_i/2) * 0.4 .* cos(phi_o - pi); 
% Eq.3
d = gamma_i - gamma_o; 
I_o = 0.5 * pi * (OD_o/2)^4 .* (1 - (g_o/OD_o)); 
I_i = 0.5 * pi * (OD_i/2)^4 .* (1 - (g_i/OD_i));

%% 2. 遍历位移 q
q_range = 5:5:40; 
colors = jet(length(q_range)); 

figure('Color', 'w', 'Name', 'CPPR Actuation Range Study');
hold on; grid on; axis equal;

for k = 1:length(q_range)
    q_target = q_range(k);
    
    % fsolve to 
    fun = @(x) solver_helper(x, q_target, n, h, d, E, I_o, I_i);
    x0 = [ones(1, n)*0.01, 1]; 
    options = optimoptions('fsolve', 'Display', 'off');
    sol = fsolve(fun, x0, options);
    
    kappa_o_sol = sol(1:n);
    
    % kinematics of CPPRs
    T_total = eye(4);
    pts = [0; 0; 0]; 
    
    for j = 1:n
        % Eq.9 & Eq.10
        kappa_j = kappa_o_sol(j) / (1 + gamma_o(j) * kappa_o_sol(j));
        l_j = h(j) / (1 - gamma_o(j) * kappa_j);
        
        theta = kappa_j * l_j;
        if abs(kappa_j) < 1e-8
            T_bj = [1,0,0,0; 0,1,0,0; 0,0,1,l_j; 0,0,0,1];
        else
            % Eq. 11
            T_bj = [cos(theta),  0,  sin(theta), (1-cos(theta))/kappa_j; ...
                    0,           1,  0,          0; ...
                   -sin(theta),  0,  cos(theta), sin(theta)/kappa_j; ...
                    0,           0,  0,          1]; 
        end
        % Eq.12
        T_rj = [1,0,0,0; 0,1,0,0; 0,0,1,c; 0,0,0,1]; % 公式 12 [cite: 282]
        % Eq.13
        T_total = T_total * T_bj * T_rj;
        pts(:, end+1) = T_total(1:3, 4);
    end
    
    % plot shape
    plot(pts(1,:), pts(3,:), 'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', ['q = ', num2str(q_target/10), ' cm']);
end

%%  plot
xlabel('X Position (mm)'); ylabel('Z Position (mm)');
title('Robot Shape Evolution (q: 0.5 to 4.0 cm)');
legend('Location', 'northeastoutside');

%% solver
function F = solver_helper(x, q_target, n, h, d, E, I_o, I_i)
    kappa_o = x(1:n);
    tau = x(n+1);
    F = zeros(1, n+1);
    for j = 1:n
        % Eq.7 静力学平衡方程 internal tension kappa is consistent across all n notches
        F(j) = (E*I_o(j)*kappa_o(j) + E*I_i(j)*(kappa_o(j)/(1 - d(j)*kappa_o(j)))) / d(j) - tau;
    end
    % Eq.8 位移约束方程 the cumulative local displacements match the global actuation input q
    F(n+1) = sum(d .* h .* kappa_o) - q_target;
end