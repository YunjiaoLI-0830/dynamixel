classdef GeometryUtils
    %GEOMETRYUTILS  Cross-section / notch geometry helpers for CPPR.
    %
    % Units: mm, MPa (consistent with the rest of this repo).
    %
    % This file contains:
    %   - Circular-segment removal math for notch modeling
    %   - Notched-tube centroid (neutral axis) offset y_bar and second moment I_x
    %   - Small utilities used by Section III-B / III-C (Fig.6)

    methods (Static)
        function [A, y_c, I_xx] = get_segment_props(r, g)
            %GET_SEGMENT_PROPS  Removed circular-segment properties.
            %
            % Inputs
            %   r : circle radius
            %   g : cut depth measured inward from the outer surface (0..2r)
            %
            % Outputs
            %   A   : area of removed segment
            %   y_c : centroid distance from circle center toward the cut
            %   I_xx: second moment of area of removed segment about the axis
            %         through the circle center parallel to the cut chord

            if g <= 0 || g >= 2*r
                A = 0; y_c = 0; I_xx = 0;
                return;
            end

            d = r - g;
            val = max(min(d/r, 1.0), -1.0);
            alpha = acos(val); % half-angle

            % Area of circular segment = sector - triangle
            area_sector   = r^2 * alpha;
            area_triangle = d * sqrt(max(r^2 - d^2, 0));
            A = area_sector - area_triangle;

            if A < 1e-12
                A = 0; y_c = 0; I_xx = 0;
                return;
            end

            % centroid of circular segment relative to circle center
            y_c = (2/3) * (r^3 * sin(alpha)^3) / A;

            % second moment of area of sector about circle center axis (parallel to chord)
            I_xx = (r^4 / 8) * (2*alpha - sin(2*alpha));
        end

        function [I_x, y_bar] = get_notched_tube_props(ro, ri, g)
            %GET_NOTCHED_TUBE_PROPS  Notched tube section properties (radii form).
            %
            % We model the remaining section as:
            %   (full annulus) - (removed outer circular segment) + (removed inner circular segment if applicable)
            %
            % Inputs
            %   ro : outer radius
            %   ri : inner radius
            %   g  : notch depth from the outer surface
            %
            % Outputs
            %   I_x   : second moment about the cut-axis, referenced to the new centroid
            %   y_bar : centroid offset relative to tube center (neutral axis offset)

            % removed segment from outer circle
            [A_cut_o, y_cut_o, I_cut_o] = GeometryUtils.get_segment_props(ro, g);

            % removed segment from inner circle (only if cut reaches inner wall)
            g_inner = g - (ro - ri);
            if g_inner <= 0
                A_cut_i = 0; y_cut_i = 0; I_cut_i = 0;
            else
                [A_cut_i, y_cut_i, I_cut_i] = GeometryUtils.get_segment_props(ri, g_inner);
            end

            % full annulus properties about circle center
            A_full_o = pi * ro^2;
            I_full_o = pi * ro^4 / 4;
            A_full_i = pi * ri^2;
            I_full_i = pi * ri^4 / 4;

            % net properties (annulus - removed segment)
            A_net = (A_full_o - A_full_i) - (A_cut_o - A_cut_i);
            if abs(A_net) < 1e-12
                I_x = 0; y_bar = 0;
                return;
            end

            % centroid (moment balance)
            moment_removed = (A_cut_o * y_cut_o) - (A_cut_i * y_cut_i);
            y_bar = -moment_removed / A_net;

            % second moment about center, then shift to centroid (parallel axis)
            I_net_center = (I_full_o - I_full_i) - (I_cut_o - I_cut_i);
            I_x = I_net_center - A_net * (y_bar^2);
        end

        function [I_x, y_bar] = get_notched_tube_props_diam(OD, ID, g)
            % wrapper: diameter form -> radii form
            [I_x, y_bar] = GeometryUtils.get_notched_tube_props(OD/2, ID/2, g);
        end

        function I_x = get_notched_I_diam(OD, ID, g)
            [I_x, ~] = GeometryUtils.get_notched_tube_props_diam(OD, ID, g);
        end

        function root = robust_fzero_scalar(f, x0, lo, hi)
            %ROBUST_FZERO_SCALAR  Scalar root finding with clamped bounds.
            opts = optimset('Display','off');

            x0 = min(max(x0, lo), hi);

            % try direct
            try
                root = fzero(f, x0, opts);
                root = min(max(root, lo), hi);
                return;
            catch
            end

            % bracket search around x0
            span = 0.1 * (hi - lo);
            for it = 1:16
                a = max(lo, x0 - span);
                b = min(hi, x0 + span);
                try
                    fa = f(a); fb = f(b);
                    if sign(fa) ~= sign(fb)
                        root = fzero(f, [a b], opts);
                        root = min(max(root, lo), hi);
                        return;
                    end
                catch
                end
                span = min(hi-lo, span * 1.8);
            end

            % fallback: minimize |f|
            g = @(x) abs(f(x));
            root = fminbnd(g, lo, hi, optimset('Display','off'));
        end

%         function gi = solve_gi_equal_EI(go, seg)
%             %SOLVE_GI_EQUAL_EI  Eq.(20): solve gi such that Eo*Io(go) = Ei*Ii(gi).
%             ODo = seg.od_o; IDo = seg.id_o;
%             ODi = seg.od_i; IDi = seg.id_i;
% 
%             Eo = seg.E; Ei = seg.E;
% 
%             Io = GeometryUtils.get_notched_I_diam(ODo, IDo, go);
%             target = (Eo/Ei) * Io;
% 
%             gi_lo = 0;
%             % maximum feasible depth: keep some backbone (avoid fully cutting through)
% %             gi_hi = max(1e-6, ODi/2 - (ODi-IDi)/2 - 1e-6);
%             gi_hi = ODi - 1e-6;
% 
% 
%             f = @(gi) GeometryUtils.get_notched_I_diam(ODi, IDi, gi) - target;
%             gi0 = min(max(go * (ODi/ODo), gi_lo), gi_hi);
% 
%             gi = GeometryUtils.robust_fzero_scalar(f, gi0, gi_lo, gi_hi);
%         end
        

          function gi = solve_gi_equal_EI(go, seg)
            %SOLVE_GI_EQUAL_EI  Eq.(20): solve gi such that Eo*Io(go) = Ei*Ii(gi).
            ODo = seg.od_o; IDo = seg.id_o;
            ODi = seg.od_i; IDi = seg.id_i;
        
            Eo = seg.E; Ei = seg.E;
        
            Io = GeometryUtils.get_notched_I_diam(ODo, IDo, go);
            target = (Eo/Ei) * Io;
        
            gi_lo = 0;
            gi_hi = ODi - 1e-6;   % g_i ∈ [0, OD_i)
        
            f = @(gi) GeometryUtils.get_notched_I_diam(ODi, IDi, gi) - target;
        
            % =====【新增：可解性检查，必须有】=====
            f_lo = f(gi_lo);
            f_hi = f(gi_hi);
            if sign(f_lo) == sign(f_hi)
                % Eq.(20) infeasible for this go:
                % even max inner cut cannot match outer stiffness
                gi = gi_lo;   % return 0 explicitly
                return;
            end
            % =====================================
        
            gi0 = min(max(go * (ODi/ODo), gi_lo), gi_hi);
            gi = GeometryUtils.robust_fzero_scalar(f, gi0, gi_lo, gi_hi);
        end


        function [alpha, d, Io, Ii, gamma_o, gamma_i] = compliance_alpha(go, gi, phi_o, seg)
            %COMPLIANCE_ALPHA  Eq.(18) compliance constant alpha.
            [Io, gamma_o, Ii, gamma_i] = seg.calculate_section_properties(go, gi, phi_o);
            d = gamma_o - gamma_i;
            alpha = d / (seg.E*Io + seg.E*Ii);
        end

        function eps_max = max_strain_notch(kappa_center, go, gi, phi_o, seg)
            %MAX_STRAIN_NOTCH  Eq.(16)-(17) maximum tensile strain in a notch pair.
            [~, gamma_o, ~, gamma_i] = seg.calculate_section_properties(go, gi, phi_o);

            ODo = seg.od_o;
            ODi = seg.od_i;
            k = abs(kappa_center);

            % Eq.(16) (outer edge a)
            eps_a_o = k * (ODo/2 - gamma_o);
            eps_a_i = k * (ODi/2 - gamma_i);

            % Eq.(17) (inner wall b formed by notch)
            eps_b_o = k * (ODo/2 + gamma_o - go);
            eps_b_i = k * (ODi/2 + gamma_i - gi);

            eps_max = max([eps_a_o, eps_b_o, eps_a_i, eps_b_i]);
        end
    end
end
