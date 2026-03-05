classdef CPPRRobot < handle
    properties
        seg     % CPPRSegment 对象
        n       % 切口数量
        h       % 切口长度数组
        g_o     % 外管切深数组
        g_i     % 内管切深数组
        phi_o   % 切口方向数组
        c       % 刚性段长度
        
        % 数值求解缓存
        last_x  % fsolve warm-start: [kappa(1:n), tau]

        % 预计算属性
        I_o
        gamma_o
        I_i
        gamma_i
        d
    end
    
    methods
        function obj = CPPRRobot(segment_params, notches)
            % 构造函数：初始化并进行几何预计算
            obj.seg = segment_params;
            obj.h = notches.h;
            obj.g_o = notches.g_o;
            obj.g_i = notches.g_i;
            obj.phi_o = notches.phi_o;
            obj.c = notches.c;
            obj.n = length(obj.h);
            
            % 初始化数组
            obj.I_o = zeros(1, obj.n);
            obj.gamma_o = zeros(1, obj.n);
            obj.I_i = zeros(1, obj.n);
            obj.gamma_i = zeros(1, obj.n);
            
            for k = 1:obj.n
                % 外管属性
                [I, y] = GeometryUtils.get_notched_tube_props(obj.seg.od_o/2, obj.seg.id_o/2, obj.g_o(k));
                obj.I_o(k) = I;
                obj.gamma_o(k) = y * cos(obj.phi_o(k));
                
                % 内管属性
                [I, y] = GeometryUtils.get_notched_tube_props(obj.seg.od_i/2, obj.seg.id_i/2, obj.g_i(k));
                obj.I_i(k) = I;
                % 内管切口方向与外管相反 (phi - pi)
                obj.gamma_i(k) = y * cos(obj.phi_o(k) - pi);
            end
            
            % Eq.3: 距离 d
            obj.d = obj.gamma_i - obj.gamma_o;
        end
        
        function F = kinematics_residual(obj, x, q_target)

            % 求解器残差函数
            % x: [kappa_o_1, ..., kappa_o_n, tau]  (n个外管曲率 + 1个内部张力)

            kappa_o = x(1:obj.n);
            tau = x(obj.n + 1);

            % 输出统一用列向量，更稳
            F = zeros(obj.n + 1, 1);

            % --- 数值保护：避免 1 - d*kappa_o 过小导致发散 ---
            denom = (1 - obj.d(:) .* kappa_o(:));
            if any(abs(denom) < 1e-8)
                % 给一个很大的残差把求解器推回可行域
                F(:) = 1e6;
                return;
            end

            % Eq.7: 力矩平衡
            term_outer = (obj.seg.E .* obj.I_o(:)) .* kappa_o(:);

            % Eq.4: 内管曲率与外管曲率关系
            kappa_i = kappa_o(:) ./ denom;
            term_inner = (obj.seg.E .* obj.I_i(:)) .* kappa_i;

            % (term_outer + term_inner)/d - tau = 0
            F(1:obj.n) = (term_outer + term_inner) ./ obj.d(:) - tau;

            % Eq.8: 几何闭环
            F(obj.n + 1) = sum(obj.d(:) .* obj.h(:) .* kappa_o(:)) - q_target;

        end
        
        function kappa_sol = solve_shape(obj, q_target)

            % 说明：外层优化会重复调用本函数，使用上一次解作为初值可显著提升稳定性与精度

            if ~isempty(obj.last_x) && numel(obj.last_x) == (obj.n + 1)
                x0 = obj.last_x(:).';  % 行向量
            else
                x0 = [ones(1, obj.n) * 1e-4, 0];
            end

            fun = @(x) obj.kinematics_residual(x, q_target);

            options = optimoptions('fsolve', ...
                'Display', 'off', ...
                'FunctionTolerance', 1e-9, ...
                'StepTolerance', 1e-9, ...
                'MaxIterations', 500, ...
                'MaxFunctionEvaluations', 5000);

            [x_sol, ~, exitflag] = fsolve(fun, x0, options);

            if exitflag <= 0
                warning('fsolve failed (exitflag=%d) for q = %.6f', exitflag, q_target);
                % 失败时不要污染缓存
            else
                obj.last_x = x_sol(:).'; % 缓存
            end

            kappa_sol = x_sol(1:obj.n);

        end 


        function obj = update_design(obj, new_h, new_g_o)

            % 更新几何参数（h 与 g_o 联合外层变量优化）
            obj.h   = new_h;
            obj.g_o = new_g_o;

            for k = 1:obj.n
                % 直接通过 Segment 计算截面参数
                [I_o, gamma_o, I_i, gamma_i] = obj.seg.calculate_section_properties( ...
                    obj.g_o(k), ...
                    obj.g_i(k), ...
                    obj.phi_o(k));

                obj.I_o(k)     = I_o;
                obj.gamma_o(k) = gamma_o;
                obj.I_i(k)     = I_i;
                obj.gamma_i(k) = gamma_i;

                obj.d(k) = gamma_i - gamma_o;
            end
        end


        function [I_o, gamma_o, I_i, gamma_i] = calculate_section_properties(obj, g_o, g_i, phi)
            % 直接委托给 CPPRSegment
            [I_o, gamma_o, I_i, gamma_i] = obj.seg.calculate_section_properties(g_o, g_i, phi);
        end

        function points = get_backbone_shape(obj, kappa_o)
            % 根据求解出的曲率重建机器人形状
            T = eye(4);
            points = zeros(3, obj.n + 1); % 预分配 (3行, n+1列)
            points(:, 1) = [0; 0; 0];
            
            for j = 1:obj.n
                % 1. 转换到中心线曲率 (Eq.9)
                kj = kappa_o(j) / (1 + obj.gamma_o(j) * kappa_o(j));
                
                % 2. 计算弯曲段弧长 (Eq.10)
                lj = obj.h(j) / (1 - obj.gamma_o(j) * kj);
                theta = kj * lj;
                
                % 3. 弯曲变换矩阵
                if abs(kj) < 1e-6 % 直线近似
                    T_b = [1, 0, 0, 0;
                           0, 1, 0, 0;
                           0, 0, 1, lj;
                           0, 0, 0, 1];
                else % 常曲率圆弧
                    s = sin(theta); c = cos(theta); inv_k = 1.0/kj;
                    % 注意：Python的变换矩阵通常是[R p]，MATLAB也是一样
                    % [ c  0  s  (1-c)/k ]
                    % [ 0  1  0     0    ]
                    % [-s  0  c   s/k    ]
                    % [ 0  0  0     1    ]
                    T_b = [c,  0,  s, inv_k*(1-c);
                           0,  1,  0,           0;
                          -s,  0,  c,     inv_k*s;
                           0,  0,  0,           1];
                end
                
                % 4. 刚性段变换矩阵 (Eq.12)
                T_r = [1, 0, 0, 0;
                       0, 1, 0, 0;
                       0, 0, 1, obj.c;
                       0, 0, 0, 1];
                   
                T = T * T_b * T_r;
                points(:, j+1) = T(1:3, 4);
            end
        end
    end
end
