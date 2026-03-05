classdef CPPRSegment < handle
    % CPPRSegment: 定义管材属性，并作为接口调用 GeometryUtils 计算截面参数
    properties
        od_o % 外管外径 [mm]
        id_o % 外管内径 [mm]
        od_i % 内管外径 [mm]
        id_i % 内管内径 [mm]
        E    % 弹性模量 [MPa]
    end
    
    methods
        % 构造函数
        function obj = CPPRSegment(od_o, id_o, od_i, id_i, E)
            obj.od_o = od_o;
            obj.id_o = id_o;
            obj.od_i = od_i;
            obj.id_i = id_i;
            obj.E = E;
        end
        
        % --- 新增：Robot 所需的接口方法 ---
        function [I_o, gamma_o, I_i, gamma_i] = calculate_section_properties(obj, g_o, g_i, phi)
            % 输入: 切深 g_o, g_i 和 切角 phi
            % 输出: 惯性矩 I 和 中性轴偏移 gamma
            
            % 1. 调用您的 GeometryUtils 计算外管属性 (半径 = 直径/2)
            [I_o_raw, y_bar_o] = GeometryUtils.get_notched_tube_props(obj.od_o/2, obj.id_o/2, g_o);
            
            % 2. 调用您的 GeometryUtils 计算内管属性
            [I_i_raw, y_bar_i] = GeometryUtils.get_notched_tube_props(obj.od_i/2, obj.id_i/2, g_i);
            
            % 3. 处理切口旋转 phi (通常平面弯曲 phi=0, cos(0)=1)
            c = cos(phi);
            gamma_o = y_bar_o * c;
            %  修改！！
            gamma_i = y_bar_i * cos(phi - pi);
%             gamma_i = y_bar_i * c;
            
            % 4. 赋值惯性矩 (简化模型：假设主轴未偏转)
            I_o = I_o_raw;
            I_i = I_i_raw;
        end
    end
end
