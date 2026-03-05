
clear; clc;

designNo = 1;  

T = readtable(fullfile('./','red_points_ordered.csv'));  % 两列 x,z

curve = struct();
curve.type = 'points';
curve.p = [T.x, T.z];   % z=-y

results = reproduce_table2_pipeline(designNo, curve);

disp(results.summary);