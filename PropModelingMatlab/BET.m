%======================================================================
% BET — Blade Element Theory de Lowry usada no TCC.
% A formulação, a discretização e as constantes do código original foram
% mantidas. Apenas os caminhos locais foram organizados.
%
% Expected optional workspace inputs:
%   geometry_file, performance_file, polar_file, target_rpm, output_dir
%
% Original script description:
%   1) Geometria (.PE0)
%   2) Bloco experimental
%   3) BET-Lowry
%   4) sMAPE acumulada (0→J) para TABELA (R², r², sMAPE_acc)
%   5) Mapas: Teórico, Referência e Teo + η(ref)
%   6) Plots comparativos η, CP, CT (teo × ref, 1 eixo)
%   7) Diagnósticos radiais (Re, Mach)
%======================================================================

clc; close all;

matlab_root = fileparts(mfilename('fullpath'));
addpath(fullfile(matlab_root, 'Functions'));

if ~exist('repo_paths','var') || isempty(repo_paths)
  repo_paths = propeller_repo_paths();
end
if ~exist('polar_file','var') || isempty(polar_file)
  polar_file = repo_paths.polar_clark_y;
end
if ~exist('output_dir','var') || isempty(output_dir)
  output_dir = repo_paths.matlab_results;
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end

%% 1) Geometria (.PE0)
run(fullfile(repo_paths.functions, 'leitura_geometria_pe0.m')); % STATION, CHORD, PITCH_LE_TE, BLADES, RADIUS ...

B = BLADES;             % número de pás
R_in = RADIUS;          % [in]
D_ft = 2*R_in/12;       % diâmetro [ft]
x_rel = STATION(:)./R_in;
n_x = numel(x_rel);

%% 2) Bloco experimental
run(fullfile(repo_paths.functions, 'leitura_desempenho_dat.m')); % J_exp, Ct_exp, Cp_exp, Pe_exp, target_rpm
J_exp = J_exp(:).';              % linha
Ct_exp = Ct_exp(:);              % coluna
Cp_exp = Cp_exp(:);
Pe_exp = Pe_exp(:);
n_J = numel(J_exp);

% aliases referência
J_ref = J_exp;
Ct_ref = Ct_exp;
Cp_ref = Cp_exp;
Eta_ref = Pe_exp;

%% 3) Banco de polares + propriedades do ar
S = load(polar_file);                 % struct com campo 'ans'
pol = S.ans.perfil;                   % NR x NM (Alpha, Cl, Cd, Cm, Re, Mach)

rho = 1.225;         % [kg/m^3]
mu = 1.81e-5;        % [Pa.s]
Vson = 343;          % [m/s]

lookup_method = 'interp'; % bilinear (Re×Mach) + linear em α

%% 4) BET-Lowry (fórmulas inalteradas)
CT_the = zeros(1, n_J);
CP_the = zeros(1, n_J);
ETA_the = zeros(1, n_J);

CL_mat = zeros(n_x, n_J);
CD_mat = zeros(n_x, n_J);
CM_mat = zeros(n_x, n_J);
Alpha_deg = zeros(n_x, n_J);
Re_mat = zeros(n_x, n_J);
Ma_mat = zeros(n_x, n_J);
Alpha_used = zeros(n_x, n_J);
Re_used = zeros(n_x, n_J);
Ma_used = zeros(n_x, n_J);

for j = 1:n_J
  J = J_exp(j);

  dCTdx = zeros(n_x, 1);
  dCPdx = zeros(n_x, 1);

  for i = 1:n_x
    c = CHORD(i);          % [in]
    r = STATION(i);        % [in]
    pitch = PITCH_LE_TE(i);% [in]

    % cinemática local
    beta = atand(pitch/(2*pi*r));     % [deg]
    phi = atand(J/(pi*x_rel(i)));     % [deg]
    alpha = beta - phi;               % [deg]

    % Re/Mach locais
    [Re_i, Ma_i] = reynolds_mach_local(r, c, J, target_rpm, R_in, rho, mu, Vson);

    % lookup (banco)
    [CL, CD, CM, a_use, Re_use, Ma_use] = lookup_clarky2(alpha, Re_i, Ma_i, pol, lookup_method);

    % diagnósticos
    Alpha_deg(i, j) = alpha;
    Re_mat(i, j) = Re_i;
    Ma_mat(i, j) = Ma_i;
    Alpha_used(i, j) = a_use;
    Re_used(i, j) = Re_use;
    Ma_used(i, j) = Ma_use;
    CL_mat(i, j) = CL;
    CD_mat(i, j) = CD;
    CM_mat(i, j) = CM;

    % Lowry
    lamT = CL*cosd(phi) - CD*sind(phi);
    lamP = CL*sind(phi) + CD*cosd(phi);

    dCTdx(i) = B*c/(4*D_ft*12).*J^2.*lamT./sind(phi).^2;
    dCPdx(i) = pi*B*c*x_rel(i)/(4*D_ft*12).*J^2.*lamP./sind(phi).^2;
  end

  CT_the(j) = trapz(x_rel, dCTdx);
  CP_the(j) = trapz(x_rel, dCPdx);
  ETA_the(j) = J*CT_the(j)/max(CP_the(j), eps); % evita inf se CP≈0
end

%% 4.1) Métricas (R², r²) + sMAPE acumulada (0→J) para a tabela
ETA_mask = ETA_the(:);
bad = ~(CT_the(:) > 0 & CP_the(:) > 0);
ETA_mask(bad) = NaN;

% métricas completas (R2, r2shape, sMAPE_mean, sMAPE_acc_curve, sMAPE_acc_final)
met_ct = series_metrics_base(Ct_exp, CT_the(:));
met_cp = series_metrics_base(Cp_exp, CP_the(:));
met_eta = series_metrics_base(Pe_exp, ETA_mask);

% curvas e valores finais de sMAPE acumulada (0→i)
smape_ct_curve = met_ct.sMAPE_acc_curve;   smape_ct_final = met_ct.sMAPE_acc_final;
smape_cp_curve = met_cp.sMAPE_acc_curve;   smape_cp_final = met_cp.sMAPE_acc_final;
smape_eta_curve = met_eta.sMAPE_acc_curve; smape_eta_final = met_eta.sMAPE_acc_final;

%% 5) Parâmetros visuais
fontName = 'Helvetica'; fs_axis = 12; fs_title = 14;
axlw = 1.0; gridA = 0.18; minorA = 0.12; figSize_cm = [14 9];

lw_ct = 2.8; lw_cp = 2.8; lw_eta = 2.8; lw_ref = 2.8;
ls_ct = '-'; ls_cp = '--'; ls_eta = '-.'; ls_ref = ':';

co_map = [0.20 0.60 0.40; 0.70 0.15 0.15; 0.10 0.25 0.55];
cp_color = [0.40 0.00 0.40];

ms_ref = 3.2; mk_ref = 'o';
cmap = @(k) co_map(mod(k-1, size(co_map, 1)) + 1, :);
lighten = @(rgb, a) (1-a)*rgb + a*[1 1 1];

format_axes = @(ax) set(ax, 'FontName', fontName, 'FontSize', fs_axis, 'Box', 'off', ...
  'TickDir', 'out', 'LineWidth', axlw, 'GridAlpha', gridA, 'MinorGridAlpha', minorA, ...
  'XMinorTick', 'on', 'YMinorTick', 'on');
newfig = @(name) set(figure('Name', name, 'Color', 'w', 'Units', 'centimeters', ...
  'Position', [1 1 figSize_cm]), 'PaperPositionMode', 'auto');

padlims = @(v, lo, hi) local_padlims(v, lo, hi);
leg = {'Location', 'best', 'Box', 'off'};

lab_CT_ref = 'C_T (ref)';
lab_CP_ref = 'C_P (ref)';
lab_ETA_ref = '\eta (ref)';

%% 6) Figuras de mapa
% Mapa Teórico
newfig('Mapa Teórico'); hold on; grid on;
yyaxis left
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cmap(2), 'DisplayName', 'C_P (teo)');
ylabel('C_T , C_P'); ylim(padlims([CT_the(:); CP_the(:); 0], 0.05, 0.08));
yyaxis right
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
ylabel('Eficiência \eta'); ylim([0 1]);
xlabel('Avanço J');
title(sprintf('Mapa Teórico – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = cmap(3); ax.YAxis(1).Color = [0.2 0.2 0.2];

% Mapa de Referência
newfig('Mapa de Referência'); hold on; grid on;
yyaxis left
plot(J_ref, Ct_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', cmap(1), 'MarkerEdgeColor', cmap(1), 'Color', lighten(cmap(1), 0.2), 'DisplayName', lab_CT_ref);
plot(J_ref, Cp_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', cp_color, 'MarkerEdgeColor', cp_color, 'Color', lighten(cp_color, 0.2), 'DisplayName', lab_CP_ref);
ylabel('C_T , C_P'); ylim(padlims([Ct_ref(:); Cp_ref(:); 0], 0.05, 0.08));
yyaxis right
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('Eficiência \eta'); ylim([0 1]);
xlabel('Avanço J');
title(sprintf('Mapa de Referência – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = [0 0 0]; ax.YAxis(1).Color = [0.2 0.2 0.2];

% Mapa Teo + η(ref)
newfig('Mapa Teo + \eta (ref)'); hold on; grid on;
yyaxis left
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cmap(2), 'DisplayName', 'C_P (teo)');
ylabel('C_T , C_P'); ylim(padlims([CT_the(:); CP_the(:); 0], 0.05, 0.08));
yyaxis right
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('Eficiência \eta'); ylim([0 1]);
xlabel('Avanço J');
title(sprintf('Mapa Teórico + \\eta (ref) – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = [0 0 0]; ax.YAxis(1).Color = [0.2 0.2 0.2];

%% 6b) Figuras comparativas simples (1 eixo)
% η × J
newfig('\eta × J – Teo × Ref'); hold on; grid on;
plot(J_exp(:), ETA_the(:), ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
plot(J_exp(:), Pe_exp(:), ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', '\eta (ref)');
ylabel('\eta'); ylim([0 1]); xlabel('Avanço J');
title(sprintf('\\eta × J – Teo × Ref  (RPM %d)', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); format_axes(gca);

% C_P × J
newfig('C_P × J – Teo × Ref'); hold on; grid on;
plot(J_exp(:), CP_the(:), ls_cp, 'LineWidth', lw_cp, 'Color', cp_color, 'DisplayName', 'C_P (teo)');
plot(J_exp(:), Cp_exp(:), ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', cp_color, 'MarkerEdgeColor', cp_color, 'Color', lighten(cp_color, 0.2), 'DisplayName', 'C_P (ref)');
ylabel('C_P'); ylim(padlims([CP_the(:); Cp_exp(:); 0], 0.05, 0.08));
xlabel('Avanço J');
title(sprintf('C_P × J – Teo × Ref  (RPM %d)', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); format_axes(gca);

% C_T × J
newfig('C_T × J – Teo × Ref'); hold on; grid on;
plot(J_exp(:), CT_the(:), ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_exp(:), Ct_exp(:), ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', cmap(1), 'MarkerEdgeColor', cmap(1), 'Color', lighten(cmap(1), 0.2), 'DisplayName', 'C_T (ref)');
ylabel('C_T'); ylim(padlims([CT_the(:); Ct_exp(:); 0], 0.05, 0.08));
xlabel('Avanço J');
title(sprintf('C_T × J – Teo × Ref  (RPM %d)', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); format_axes(gca);

%% 7) Tabela resumo (console + CSV opcional)
labels = {'C_T', 'C_P', '\eta'};
csvname = fullfile(output_dir, sprintf('metrics_bet_RPM_%d.csv', target_rpm));
print_metrics_table_acc(met_ct, met_cp, met_eta, ...
  [smape_ct_final, smape_cp_final, smape_eta_final], labels, target_rpm, csvname);

%% 8) Diagnóstico radial (Re/Mach)
co_rad = [0.20 0.60 0.40; 0.20 0.40 0.70; 0.50 0.40 0.70; 0.65 0.45 0.25; 0.90 0.70 0.30];
crad = @(k) co_rad(mod(k-1, size(co_rad, 1)) + 1, :);

nJshow = min(6, n_J);  nXshow = min(6, n_x);
idxJ = unique(round(linspace(1, n_J, nJshow)));
idxX = unique(round(linspace(1, n_x, nXshow)));
legJ = arrayfun(@(jj) sprintf('J=%.3g', J_exp(jj)), idxJ, 'UniformOutput', false);
legX = arrayfun(@(ii) sprintf('x=%.2f', x_rel(ii)), idxX, 'UniformOutput', false);

% Re(x)
newfig('Re(x)'); hold on; grid on;
for k = 1:numel(idxJ), j = idxJ(k);
  semilogy(x_rel, Re_used(:, j), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('x = r/R'); ylabel('Re'); title('Reynolds ao longo da pá');
legend(legJ, 'Location', 'best', 'Box', 'off'); format_axes(gca);

% Re(J)
newfig('Re(J)'); hold on; grid on;
for k = 1:numel(idxX), i = idxX(k);
  semilogy(J_exp(:), Re_used(i, :).', '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('J'); ylabel('Re'); title('Reynolds em estações fixas');
legend(legX, 'Location', 'best', 'Box', 'off'); format_axes(gca);

% Mach(x)
newfig('Mach(x)'); hold on; grid on;
for k = 1:numel(idxJ), j = idxJ(k);
  plot(x_rel, Ma_used(:, j), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('x = r/R'); ylabel('Mach'); title('Mach ao longo da pá');
legend(legJ, 'Location', 'best', 'Box', 'off'); format_axes(gca);

% Mach(J)
newfig('Mach(J)'); hold on; grid on;
for k = 1:numel(idxX), i = idxX(k);
  plot(J_exp(:), Ma_used(i, :).', '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('J'); ylabel('Mach'); title('Mach em estações fixas');
legend(legX, 'Location', 'best', 'Box', 'off'); format_axes(gca);
