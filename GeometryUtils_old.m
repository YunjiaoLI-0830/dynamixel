classdef GeometryUtils_old
    methods (Static)
        function [A, y_c, I_xx] = get_segment_props(r, g)
            % 计算圆弓形（被切掉的部分）的截面属性
            % input: r=半径, g=切深
            % output: A=面积, y_c=形心距离圆心的距离, I_xx=关于平行于切线轴的惯性矩
            
            if g >= 2 * r % 全切
                A = 0; y_c = 0; I_xx = 0; return;
            end
            if g <= 0     % 没切
                A = 0; y_c = 0; I_xx = 0; return;
            end
            
            d = r - g;
            % alpha 是扇形半角
            % 强制将结果保存在[-1,1]之间，避免运算报错
            val = max(min(d/r, 1.0), -1.0);
            alpha = acos(val);
            
            % 弓形面积 = 扇形 - 三角形
            area_sector = (r^2) * alpha;
            area_triangle = d * sqrt(r^2 - d^2);
            A = area_sector - area_triangle;
            
            if A < 1e-9
                A = 0; y_c = 0; I_xx = 0; return;
            end
            
            % 弓形形心 (相对于圆心)
            y_c = (2/3) * (r^3 * sin(alpha)^3) / A;
            
            % 弓形惯性矩 (关于圆心且平行于切线的轴)
            % 公式: (R^4 / 8) * (2*alpha - sin(2*alpha))
            I_xx = (r^4 / 8) * (2*alpha - sin(2*alpha));
        end
        
        function [I_x, y_bar] = get_notched_tube_props(ro, ri, g)
            % 核心逻辑：剩余截面 = (完整外圆 - 完整内圆) - (外切弓形 - 内切弓形)
            % input: r_o外径，r_i内径，g=切深
            % output: I_x绕新的形心的惯性矩, y_bar新的形心
            
            % 1. 计算被切掉的“外弓形”属性
            [A_cut_o, y_cut_o, I_cut_o] = GeometryUtils.get_segment_props(ro, g);
            
            % 2. 计算被切掉的“内弓形”属性
            g_inner = g - (ro - ri);
            if g_inner <= 0
                A_cut_i = 0; y_cut_i = 0; I_cut_i = 0;
            else
                [A_cut_i, y_cut_i, I_cut_i] = GeometryUtils.get_segment_props(ri, g_inner);
            end
            
            % 3. 完整管材属性
            A_full_o = pi * ro^2;
            I_full_o = pi * ro^4 / 4;
            A_full_i = pi * ri^2;
            I_full_i = pi * ri^4 / 4;
            
            % 4. 计算剩余部分的净属性 (Net Properties)
            % 净面积
            A_net = (A_full_o - A_full_i) - (A_cut_o - A_cut_i);
            
            % 净形心 (Neutral Axis Offset gamma)
            % Moment_net = 0 - (M_cut_o - M_cut_i)
            moment_removed = (A_cut_o * y_cut_o) - (A_cut_i * y_cut_i);
            y_bar = -moment_removed / A_net;
            
            % 净惯性矩 (关于圆心)
            I_net_center = (I_full_o - I_full_i) - (I_cut_o - I_cut_i);
            
            % 平行轴定理移轴到新形心: I_new = I_center - A * d^2
            % 需要计算绕新的形心的惯性矩，不能算绕圆心的惯性矩
            I_x = I_net_center - A_net * (y_bar^2);
        end
    end
end