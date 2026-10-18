%======================================================================
% BEMT — Blade Element Momentum Theory usada no TCC.
% A formulação, a discretização e as constantes do código original foram
% mantidas. Apenas os caminhos locais foram organizados.
%
% Reproducibility note: masterEq.m retains the Lowry analytical polar used
% by the original angle solver; final loads use the Clark-Y database.
%
% Expected optional workspace inputs:
%   geometry_file, performance_file, polar_file, target_rpm, output_dir
%
% Original script description:
%   1) Lê geometria e dados de referência
%   2) Calcula mapa teórico com elemento de pá + momentum 
%   3) Métricas globais (R², r² e sMAPE)
%   4) Mapas (Teórico, Referência, Teo+ηref) e comparativos
%   5) Diagnósticos radiais (Re, Mach)
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

%% 1) GEOMETRIA
run(fullfile(repo_paths.functions, 'leitura_geometria_pe0.m')); % STATION, CHORD, PITCH_LE_TE, RADIUS, BLADES

B = BLADES;
R_in = RADIUS;                 % [in]
D_in = 2*R_in;                 % [in]

x_rel = STATION(:)./R_in;      % r/R (hub já excluído no .PE0)
n_x = numel(x_rel);

%% 2) DADOS DE REFERÊNCIA (FABRICANTE)
run(fullfile(repo_paths.functions, 'leitura_desempenho_dat.m')); % J_exp, Ct_exp, Cp_exp, Pe_exp, target_rpm

J_ref = J_exp(:)';      % mantém pontos de J do experimental
Ct_ref = Ct_exp(:);
Cp_ref = Cp_exp(:);
Eta_ref = Pe_exp(:);
n_J = numel(J_ref);

%% 3) PREPARO (POLARES + AR)
S = load(polar_file);                 % grid de structs em S.ans.perfil
pol = S.ans.perfil;

rho = 1.225;   % kg/m^3
mu = 1.81e-5;  % Pa·s
Vson = 343;    % m/s

lookup_method = 'interp';

CT_the = zeros(n_J, 1);
CP_the = zeros(n_J, 1);
ETA_the = zeros(n_J, 1);

deg2rad = pi/180;
optFzero = optimset('Display', 'off');

alpha_grid = -10:0.5:25;  % grade de busca do solver

% diagnósticos
CL_x = NaN(n_x, n_J);
CD_x = NaN(n_x, n_J);
ALF_x = NaN(n_x, n_J);
RE_x = NaN(n_x, n_J);
MA_x = NaN(n_x, n_J);

%% 4) LOOP SOBRE J  [núcleo inalterado]
for j = 1:n_J
  Jt = J_ref(j);

  % ignora apenas o teórico para J < 0.05 (mantém dados de referência)
  if Jt < 0.05
    CT_the(j) = NaN; CP_the(j) = NaN; ETA_the(j) = NaN;
    CL_x(:, j) = NaN; CD_x(:, j) = NaN; ALF_x(:, j) = NaN;
    RE_x(:, j) = NaN; MA_x(:, j) = NaN;
    continue
  end

  dCTdx = zeros(n_x, 1);
  dCPdx = zeros(n_x, 1);

  for i = 1:n_x
    x = x_rel(i);
    c = CHORD(i);
    r = STATION(i);
    beta = atand(PITCH_LE_TE(i)/(2*pi*r));
    sigma = B*c/(2*pi*x*D_in/2);

    % resolve α
    fvals = arrayfun(@(a) masterEq(a, beta, Jt, x, sigma, B), alpha_grid);
    idx = find(diff(sign(fvals)) ~= 0, 1);
    alpha = NaN;
    if ~isempty(idx)
      a_lo = alpha_grid(idx);
      a_hi = alpha_grid(idx+1);
      try
        alpha = fzero(@(a) masterEq(a, beta, Jt, x, sigma, B), [a_lo a_hi], optFzero);
      catch
        alpha = NaN;
      end
    end

    if isnan(alpha)
      dCTdx(i) = 0; dCPdx(i) = 0;
      continue
    end

    phi = beta - alpha;    phi_r = phi*deg2rad;
    if abs(sin(phi_r)) < 1e-9
      dCTdx(i) = 0; dCPdx(i) = 0;
      continue
    end

    % Re/Mach locais + polares
    [Re_i, Ma_i] = reynolds_mach_local(r, c, Jt, target_rpm, R_in, rho, mu, Vson);
    [CL, CD] = lookup_clarky2(alpha, Re_i, Ma_i, pol, lookup_method);

    % diagnósticos
    CL_x(i, j) = CL;
    CD_x(i, j) = CD;
    ALF_x(i, j) = alpha;
    RE_x(i, j) = Re_i;
    MA_x(i, j) = Ma_i;

    % Lowry-Momentum
    lamT = CL*cosd(phi) - CD*sind(phi);
    lamP = CL*sind(phi) + CD*cosd(phi);

    F = prandtlF(B, x, phi_r);
    a_ax = 1/(4*F*sin(phi_r).^2/(lamT*sigma) - 1);
    a_tan = 1/(4*F*sin(phi_r).*cos(phi_r)/(lamP*sigma) + 1); %#ok<NASGU>

    fact = (1 + a_ax).^2.*(1./sind(phi)).^2;
    dCTdx(i) = 0.25*Jt^2*fact*c*B*lamT/D_in;
    dCPdx(i) = 0.25*pi*x*Jt^2*fact*c*B*lamP/D_in;
  end

  CT_the(j) = trapz(x_rel, dCTdx);
  CP_the(j) = trapz(x_rel, dCPdx);
  ETA_the(j) = Jt*CT_the(j)/max(CP_the(j), eps);
end

%% 5) MÉTRICAS GLOBAIS (R², r² e sMAPE) — TABELA
ETA_mask = ETA_the;
bad = ~(CT_the > 0 & CP_the > 0);
ETA_mask(bad) = NaN;

% usa a base unificada: R2, r2shape, sMAPE_mean, sMAPE_acc_curve, sMAPE_acc_final
met_ct  = series_metrics_base(Ct_ref , CT_the );
met_cp  = series_metrics_base(Cp_ref , CP_the );
met_eta = series_metrics_base(Eta_ref, ETA_mask);


%% 6) PLOTS — parâmetros visuais
fontName = 'Helvetica'; fs_axis = 12; fs_title = 14;
axlw = 1.0; gridA = 0.18; minorA = 0.12; figSize_cm = [14 9];

lw_ct = 2.8; lw_cp = 2.8; lw_eta = 2.8; lw_ref = 2.8;
ms_ref = 3.2; mk_ref = 'o';

ls_ct = '-'; ls_cp = '--'; ls_eta = '-.'; ls_ref = ':';

co_map = [0.20 0.60 0.40; 0.70 0.15 0.15; 0.10 0.25 0.55];
cp_color = [0.40 0.00 0.40];

cmap = @(k) co_map(mod(k-1, size(co_map, 1)) + 1, :);
lighten = @(rgb, a) (1 - a)*rgb + a*[1 1 1];

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

%% FIGURAS DE MAPA
% Mapa Teórico
newfig('Mapa Teórico'); hold on; grid on;
yyaxis left
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cmap(2), 'DisplayName', 'C_P (teo)');
ylabel('C_T , C_P'); ylim(padlims([CT_the; CP_the; 0], 0.05, 0.08));
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
ylabel('C_T , C_P'); ylim(padlims([Ct_ref; Cp_ref; 0], 0.05, 0.08));
yyaxis right
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('Eficiência \eta'); ylim([0 1]);
xlabel('Avanço J');
title(sprintf('Mapa de Referência – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = [0 0 0]; ax.YAxis(1).Color = [0.2 0.2 0.2];

% Mapa Teo + η (ref)
newfig('Mapa Teo + \eta (ref)'); hold on; grid on;
yyaxis left
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cmap(2), 'DisplayName', 'C_P (teo)');
ylabel('C_T , C_P'); ylim(padlims([CT_the; CP_the; 0], 0.05, 0.08));
yyaxis right
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('Eficiência \eta'); ylim([0 1]);
xlabel('Avanço J');
title(sprintf('Mapa Teórico + \\eta (ref) – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName);
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = [0 0 0]; ax.YAxis(1).Color = [0.2 0.2 0.2];

%% COMPARATIVOS (1 eixo)
% C_T × J
newfig('C_T × J – Teo × Ref'); hold on; grid on;
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, Ct_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', cmap(1), 'MarkerEdgeColor', cmap(1), 'Color', lighten(cmap(1), 0.2), 'DisplayName', lab_CT_ref);
ylabel('C_T'); ylim(padlims([CT_the; Ct_ref; 0], 0.05, 0.08));
xlabel('Avanço J');
title(sprintf('C_T × J – Teo × Ref  (RPM %d)', target_rpm));
legend(leg{:}); format_axes(gca);

% C_P × J
newfig('C_P × J – Teo × Ref'); hold on; grid on;
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cp_color, 'DisplayName', 'C_P (teo)');
plot(J_ref, Cp_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', cp_color, 'MarkerEdgeColor', cp_color, 'Color', lighten(cp_color, 0.2), 'DisplayName', lab_CP_ref);
ylabel('C_P'); ylim(padlims([CP_the; Cp_ref; 0], 0.05, 0.08));
xlabel('Avanço J');
title(sprintf('C_P × J – Teo × Ref  (RPM %d)', target_rpm));
legend(leg{:}); format_axes(gca);

% η × J
newfig('\eta × J – Teo × Ref'); hold on; grid on;
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('\eta'); ylim([0 1]);
xlabel('Avanço J');
title(sprintf('\\eta × J – Teo × Ref  (RPM %d)', target_rpm));
legend(leg{:}); format_axes(gca);

%% TABELA (console + CSV opcional)
labels = {'C_T','C_P','\eta'};
csvname = fullfile(output_dir, sprintf('metrics_bemt_RPM_%d.csv', target_rpm));

% usa os campos da base unificada
smape_ct_final  = met_ct.sMAPE_acc_final;
smape_cp_final  = met_cp.sMAPE_acc_final;
smape_eta_final = met_eta.sMAPE_acc_final;

print_metrics_table_acc( ...
    met_ct, met_cp, met_eta, ...
    [smape_ct_final smape_cp_final smape_eta_final], ...
    labels, target_rpm, csvname);

%% DIAGNÓSTICOS RADIAIS (Re/Mach)
co_rad = [0.20 0.60 0.40; 0.20 0.40 0.70; 0.50 0.40 0.70; 0.65 0.45 0.25; 0.90 0.70 0.30];
crad = @(k) co_rad(mod(k-1, size(co_rad, 1)) + 1, :);

nJshow = min(6, n_J);  nXshow = min(6, n_x);
idxJ = unique(round(linspace(1, n_J, nJshow)));
idxX = unique(round(linspace(1, n_x, nXshow)));
legJ = arrayfun(@(jj) sprintf('J=%.3g', J_ref(jj)), idxJ, 'UniformOutput', false);
legX = arrayfun(@(ii) sprintf('x=%.2f', x_rel(ii)), idxX, 'UniformOutput', false);

% Re(x)
newfig('Re(x)'); hold on; grid on;
for k = 1:numel(idxJ), j = idxJ(k)
  semilogy(x_rel, RE_x(:, j), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('x = r/R'); ylabel('Re'); title('Reynolds ao longo da pá');
legend(legJ, 'Location', 'best', 'Box', 'off'); format_axes(gca);

% Re(J)
newfig('Re(J)'); hold on; grid on;
for k = 1:numel(idxX), i = idxX(k)
  semilogy(J_ref, RE_x(i, :), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('J'); ylabel('Re'); title('Reynolds em estações fixas');
legend(legX, 'Location', 'best', 'Box', 'off'); format_axes(gca);

% Mach(x)
newfig('Mach(x)'); hold on; grid on;
for k = 1:numel(idxJ), j = idxJ(k)
  plot(x_rel, MA_x(:, j), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('x = r/R'); ylabel('Mach'); title('Mach ao longo da pá');
legend(legJ, 'Location', 'best', 'Box', 'off'); format_axes(gca);

% Mach(J)
newfig('Mach(J)'); hold on; grid on;
for k = 1:numel(idxX), i = idxX(k)
  plot(J_ref, MA_x(i, :), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('J'); ylabel('Mach'); title('Mach em estações fixas');
legend(legX, 'Location', 'best', 'Box', 'off'); format_axes(gca);
