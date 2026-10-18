%COMPAREMETHODS Compara BET, BEMT e Vortex no mesmo caso APC 22x10.
% Este arquivo apenas executa os três modelos; não altera suas equações.
clear; clc; close all;

matlab_root = fileparts(mfilename('fullpath'));
addpath(fullfile(matlab_root, 'Functions'));
repo_paths = propeller_repo_paths();

geometry_file = repo_paths.geometry_22x10;
performance_file = repo_paths.performance_22x10;
polar_file = repo_paths.polar_clark_y;
target_rpm = 5000;
plot_reference_data = false;
output_dir = repo_paths.matlab_results;

run(fullfile(matlab_root, 'BET.m'));
comparison.reference = struct('J', J_ref(:), 'CT', Ct_ref(:), ...
    'CP', Cp_ref(:), 'eta', Eta_ref(:));
comparison.BET = struct('J', J_ref(:), 'CT', CT_the(:), ...
    'CP', CP_the(:), 'eta', ETA_the(:));

run(fullfile(matlab_root, 'BEMT.m'));
comparison.BEMT = struct('J', J_ref(:), 'CT', CT_the(:), ...
    'CP', CP_the(:), 'eta', ETA_the(:));

clear busca_polares_corretas;
run(fullfile(matlab_root, 'Vortex.m'));
comparison.Vortex = struct('J', J_ref(:), 'CT', CT_the(:), ...
    'CP', CP_the(:), 'eta', ETA_the(:));

colors = struct('BET', [0.35 0.38 0.42], ...
    'BEMT', [0.12 0.48 0.66], 'Vortex', [0.86 0.44 0.08]);
figure('Color', 'w', 'Name', 'BET-BEMT-Vortex comparison');
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
fields = {'CT', 'CP', 'eta'};
labels = {'C_T', 'C_P', '\eta'};
for k = 1:numel(fields)
    nexttile; hold on; grid on;
    name = fields{k};
    plot(comparison.reference.J, comparison.reference.(name), 'k:o', ...
        'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'APC reference');
    plot(comparison.BET.J, comparison.BET.(name), '-', ...
        'Color', colors.BET, 'LineWidth', 1.8, 'DisplayName', 'BET');
    plot(comparison.BEMT.J, comparison.BEMT.(name), '--', ...
        'Color', colors.BEMT, 'LineWidth', 1.8, 'DisplayName', 'BEMT');
    plot(comparison.Vortex.J, comparison.Vortex.(name), '-.', ...
        'Color', colors.Vortex, 'LineWidth', 1.8, 'DisplayName', 'Vortex');
    xlabel('Advance ratio, J'); ylabel(labels{k});
    if strcmp(name, 'eta'), ylim([0 1]); end
end
legend('Location', 'southoutside', 'Orientation', 'horizontal');
title(t, sprintf('APC 22x10 at %d rpm', target_rpm));
