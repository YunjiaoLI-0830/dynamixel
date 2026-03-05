function results = reproduce_table2_pipeline(designNo, curve)
%REPRODUCE_TABLE2_PIPELINE  Run Fig.6 design using Table II parameters.
%
% Usage:
%   curve.type = 'segments'; curve.segments = [...];
%   results = reproduce_table2_pipeline(1, curve);
%
% Inputs
%   designNo : integer 1..6
%   curve    : desired curve struct (see sample_curve_for_fig6)
%
% Output results struct:
%   .params   : design params used (Table II + derived tube IDs)
%   .out      : output from design_fig6_from_curve
%   .summary  : min/max of h and g_o similar to Table II

if nargin < 2
    error('Usage: results = reproduce_table2_pipeline(designNo, curve)');
end
if ~isscalar(designNo) || designNo < 1 || designNo > 6
    error('designNo must be an integer 1..6.');
end

T = table2_params();
row = T(designNo);

% Convert Table II values to model units
ODo = row.ODo_mm;
ODi = row.ODi_mm;
to = 1; ti = 1; % Table II note: to = ti = 1 mm

IDo = ODo - 2*to;
IDi = ODi - 2*ti;

E_MPa = row.E_GPa * 1000;
epsFS = row.epsFS_percent / 100; % fraction

seg = CPPRSegment(ODo, IDo, ODi, IDi, E_MPa);

params = struct();
params.Lcurve = row.Lcurve_mm;
params.beta   = row.beta;
params.n      = row.n;
params.qmax   = row.qmax_mm
params.epsFS  = epsFS;
params.segment_params = seg;
params.g_bounds = [0.2, params.segment_params.od_o - 0.2];

params.method = "sequential";

out = design_fig6_from_curve(curve, params);

results = struct();
results.params = params;
results.out = out;

results.summary = struct();
results.summary.h_min = min(out.notches.h);
results.summary.h_max = max(out.notches.h);
results.summary.go_min = min(out.notches.g_o);
results.summary.go_max = max(out.notches.g_o);
results.summary.max_strain = out.opt.max_strain;
results.summary.epsFS = epsFS;

end

function T = table2_params()
% Values transcribed from Table II (page 10 of the paper).
% Columns: E(GPa), eps_max(%), eps_FS(%), ODo(mm), ODi(mm), Lcurve(mm), beta, n, qmax(mm)
T = struct([]);

T(1).E_GPa = 1.58;  T(1).epsmax_percent = 7.00;  T(1).epsFS_percent = 5.80;  T(1).ODo_mm = 8.0;  T(1).ODi_mm = 5.2;  T(1).Lcurve_mm = 120; T(1).beta = 0.5;  T(1).n = 10; T(1).qmax_mm = 13.92;
T(2).E_GPa = 1.58;  T(2).epsmax_percent = 7.00;  T(2).epsFS_percent = 6.36;  T(2).ODo_mm = 6.0;  T(2).ODi_mm = 2.8;  T(2).Lcurve_mm = 100; T(2).beta = 0.55; T(2).n = 12; T(2).qmax_mm = 13.44;
T(3).E_GPa = 1.80;  T(3).epsmax_percent = 4.50;  T(3).epsFS_percent = 3.00;  T(3).ODo_mm = 10.0; T(3).ODi_mm = 6.54; T(3).Lcurve_mm = 150; T(3).beta = 0.5;  T(3).n = 10; T(3).qmax_mm = 16.26;
T(4).E_GPa = 0.502; T(4).epsmax_percent = 18.0;  T(4).epsFS_percent = 6.67;  T(4).ODo_mm = 8.0;  T(4).ODi_mm = 4.8;  T(4).Lcurve_mm = 120; T(4).beta = 0.5;  T(4).n = 12; T(4).qmax_mm = 13.22;
T(5).E_GPa = 0.502; T(5).epsmax_percent = 18.0;  T(5).epsFS_percent = 6.67;  T(5).ODo_mm = 8.0;  T(5).ODi_mm = 4.3;  T(5).Lcurve_mm = 150; T(5).beta = 0.65; T(5).n = 15; T(5).qmax_mm = 17.96;
T(6).E_GPa = 1.80;  T(6).epsmax_percent = 4.50;  T(6).epsFS_percent = 3.91;  T(6).ODo_mm = 8.0;  T(6).ODi_mm = 4.8;  T(6).Lcurve_mm = 150; T(6).beta = 0.45; T(6).n = 15; T(6).qmax_mm = 15.75;

end
