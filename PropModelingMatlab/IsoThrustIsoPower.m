%% ============================================================
% Compara a hélice MEL e a APC 22x10 em isopotência e isotração.
% Os valores e a interpolação do código do TCC foram mantidos.
%
% APC 22x10 @ 4300 rpm REAL
% Comparacoes:
%   1) Mesmo J
%   2) Iso-potencia
%   3) Iso-tracao
%
% Saida completa no Command Window para montar tabela do artigo
%% ============================================================

clearvars -except DAT_file
clc; close all;

matlab_root = fileparts(mfilename('fullpath'));
addpath(fullfile(matlab_root, 'Functions'));

%% ===================== DADOS FIXOS ==========================
rho = 1.225;              % [kg/m^3]
D   = 0.559;              % [m] usar o MESMO do seu codigo principal
RPM = 4300;               % [rpm]
n   = RPM/60;             % [1/s]

% ponto da sua helice otimizada
J_ref    = 0.297;
CT_ref   = 0.0464;
CP_ref   = 0.020352;      % equivalente a 500 W em 4300 rpm
eta_ref  = 0.6780;
T_ref    = 28.49;         % [N]
P_ref    = 500.0;         % [W]
Cl_ref   = 0.760;
zeta_ref = 0.6807;
AF_ref   = 104.5;

fprintf('================ DADOS DE REFERENCIA ================\n');
fprintf('Helice otimizada:\n');
fprintf('Cl*      = %.4f\n', Cl_ref);
fprintf('zeta     = %.4f\n', zeta_ref);
fprintf('J        = %.4f\n', J_ref);
fprintf('CT       = %.6f\n', CT_ref);
fprintf('CP       = %.6f\n', CP_ref);
fprintf('eta      = %.4f\n', eta_ref);
fprintf('T        = %.3f N\n', T_ref);
fprintf('P        = %.3f W\n', P_ref);
fprintf('AF       = %.1f\n\n', AF_ref);

%% ===================== LER DADOS APC ========================
if ~exist('DAT_file','var') || isempty(DAT_file)
    repo_paths = propeller_repo_paths();
    DAT_file = repo_paths.performance_22x10;
end

EXP4 = read_APC_PER3_nearest(DAT_file, 4000);
EXP5 = read_APC_PER3_nearest(DAT_file, 5000);

J4  = EXP4.J(:);   Ct4 = EXP4.Ct(:);   Cp4 = EXP4.Cp(:);   Pe4 = EXP4.Pe(:);
J5  = EXP5.J(:);   Ct5 = EXP5.Ct(:);   Cp5 = EXP5.Cp(:);   Pe5 = EXP5.Pe(:);

fprintf('Curvas APC usadas: 4000 rpm e 5000 rpm\n');

%% ===================== MALHA COMUM DE J =====================
Jmin = max(min(J4), min(J5));
Jmax = min(max(J4), max(J5));
J_common = linspace(Jmin, Jmax, 1000)';

%% ===================== INTERPOLAR EM J ======================
Ct4i = interp1(J4, Ct4, J_common, 'pchip');
Cp4i = interp1(J4, Cp4, J_common, 'pchip');
Pe4i = interp1(J4, Pe4, J_common, 'pchip');

Ct5i = interp1(J5, Ct5, J_common, 'pchip');
Cp5i = interp1(J5, Cp5, J_common, 'pchip');
Pe5i = interp1(J5, Pe5, J_common, 'pchip');

%% ===================== INTERPOLAR PARA 4300 RPM ============
w = (RPM - 4000)/(5000 - 4000);   % 0.3

Ct4300 = Ct4i + w*(Ct5i - Ct4i);
Cp4300 = Cp4i + w*(Cp5i - Cp4i);
Pe4300 = Pe4i + w*(Pe5i - Pe4i);

%% ===================== CURVAS EM UNIDADES FISICAS ==========
T4300 = Ct4300 .* rho .* n.^2 .* D.^4;   % [N]
P4300 = Cp4300 .* rho .* n.^3 .* D.^5;   % [W]
V4300 = J_common .* n .* D;              % [m/s]

%% ===================== 1) MESMO J ===========================
J_same = J_ref;

CT_same  = interp1(J_common, Ct4300, J_same, 'pchip');
CP_same  = interp1(J_common, Cp4300, J_same, 'pchip');
eta_same = interp1(J_common, Pe4300, J_same, 'pchip');
T_same   = interp1(J_common, T4300,  J_same, 'pchip');
P_same   = interp1(J_common, P4300,  J_same, 'pchip');
V_same   = interp1(J_common, V4300,  J_same, 'pchip');

%% ===================== 2) ISO-POTENCIA ======================
% Procurar CP = CP_ref no ramo descendente
[~, idx_peak_cp] = max(Cp4300);

J_branch_cp  = J_common(idx_peak_cp:end);
Cp_branch    = Cp4300(idx_peak_cp:end);
Ct_branch_cp = Ct4300(idx_peak_cp:end);
Pe_branch_cp = Pe4300(idx_peak_cp:end);
T_branch_cp  = T4300(idx_peak_cp:end);
P_branch_cp  = P4300(idx_peak_cp:end);
V_branch_cp  = V4300(idx_peak_cp:end);

[Cp_sorted, idx_cp] = sort(Cp_branch);
J_sorted_cp   = J_branch_cp(idx_cp);
Ct_sorted_cp  = Ct_branch_cp(idx_cp);
Pe_sorted_cp  = Pe_branch_cp(idx_cp);
T_sorted_cp   = T_branch_cp(idx_cp);
P_sorted_cp   = P_branch_cp(idx_cp);
V_sorted_cp   = V_branch_cp(idx_cp);

[Cp_unique, ia_cp] = unique(Cp_sorted);
J_unique_cp   = J_sorted_cp(ia_cp);
Ct_unique_cp  = Ct_sorted_cp(ia_cp);
Pe_unique_cp  = Pe_sorted_cp(ia_cp);
T_unique_cp   = T_sorted_cp(ia_cp);
P_unique_cp   = P_sorted_cp(ia_cp);
V_unique_cp   = V_sorted_cp(ia_cp);

J_isoP   = interp1(Cp_unique, J_unique_cp,  CP_ref, 'pchip');
CT_isoP  = interp1(Cp_unique, Ct_unique_cp, CP_ref, 'pchip');
CP_isoP  = interp1(Cp_unique, Cp_unique,    CP_ref, 'pchip');
eta_isoP = interp1(Cp_unique, Pe_unique_cp, CP_ref, 'pchip');
T_isoP   = interp1(Cp_unique, T_unique_cp,  CP_ref, 'pchip');
P_isoP   = interp1(Cp_unique, P_unique_cp,  CP_ref, 'pchip');
V_isoP   = interp1(Cp_unique, V_unique_cp,  CP_ref, 'pchip');

%% ===================== 3) ISO-TRACAO ========================
[T_sorted, idx_t] = sort(T4300);
J_sorted_t   = J_common(idx_t);
Ct_sorted_t  = Ct4300(idx_t);
Cp_sorted_t  = Cp4300(idx_t);
Pe_sorted_t  = Pe4300(idx_t);
P_sorted_t   = P4300(idx_t);
V_sorted_t   = V4300(idx_t);

[T_unique, ia_t] = unique(T_sorted);
J_unique_t   = J_sorted_t(ia_t);
Ct_unique_t  = Ct_sorted_t(ia_t);
Cp_unique_t  = Cp_sorted_t(ia_t);
Pe_unique_t  = Pe_sorted_t(ia_t);
P_unique_t   = P_sorted_t(ia_t);
V_unique_t   = V_sorted_t(ia_t);

J_isoT   = interp1(T_unique, J_unique_t,  T_ref, 'pchip');
CT_isoT  = interp1(T_unique, Ct_unique_t, T_ref, 'pchip');
CP_isoT  = interp1(T_unique, Cp_unique_t, T_ref, 'pchip');
eta_isoT = interp1(T_unique, Pe_unique_t, T_ref, 'pchip');
P_isoT   = interp1(T_unique, P_unique_t,  T_ref, 'pchip');
V_isoT   = interp1(T_unique, V_unique_t,  T_ref, 'pchip');

%% ===================== AF DA APC ============================
% manter o valor que voce ja obteve no codigo principal
AF_apc = 120.9;

%% ===================== PRINT FINAL ==========================
fprintf('\n================ RESULTADOS FINAIS APC @ 4300 rpm ================\n');

fprintf('\n--- 1) APC NO MESMO J DA HELICE OTIMIZADA ---\n');
fprintf('J        = %.6f\n', J_same);
fprintf('V [m/s]  = %.3f\n', V_same);
fprintf('CT       = %.6f\n', CT_same);
fprintf('CP       = %.6f\n', CP_same);
fprintf('eta      = %.4f\n', eta_same);
fprintf('T [N]    = %.3f\n', T_same);
fprintf('P [W]    = %.3f\n', P_same);
fprintf('AF       = %.1f\n', AF_apc);

fprintf('\n--- 2) APC EM ISO-POTENCIA (CP = %.6f / P = %.1f W) ---\n', CP_ref, P_ref);
fprintf('J        = %.6f\n', J_isoP);
fprintf('V [m/s]  = %.3f\n', V_isoP);
fprintf('CT       = %.6f\n', CT_isoP);
fprintf('CP       = %.6f\n', CP_isoP);
fprintf('eta      = %.4f\n', eta_isoP);
fprintf('T [N]    = %.3f\n', T_isoP);
fprintf('P [W]    = %.3f\n', P_isoP);
fprintf('AF       = %.1f\n', AF_apc);

fprintf('\n--- 3) APC EM ISO-TRACAO (T = %.3f N) ---\n', T_ref);
fprintf('J        = %.6f\n', J_isoT);
fprintf('V [m/s]  = %.3f\n', V_isoT);
fprintf('CT       = %.6f\n', CT_isoT);
fprintf('CP       = %.6f\n', CP_isoT);
fprintf('eta      = %.4f\n', eta_isoT);
fprintf('T [N]    = %.3f\n', T_ref);
fprintf('P [W]    = %.3f\n', P_isoT);
fprintf('AF       = %.1f\n', AF_apc);

fprintf('\n================ COMPARACAO COM HELICE OTIMIZADA ================\n');
fprintf('Helice otimizada : J=%.6f | CT=%.6f | CP=%.6f | eta=%.4f | T=%.3f N | P=%.3f W | AF=%.1f\n', ...
    J_ref, CT_ref, CP_ref, eta_ref, T_ref, P_ref, AF_ref);

fprintf('APC mesmo J      : J=%.6f | CT=%.6f | CP=%.6f | eta=%.4f | T=%.3f N | P=%.3f W | AF=%.1f\n', ...
    J_same, CT_same, CP_same, eta_same, T_same, P_same, AF_apc);

fprintf('APC iso-potencia : J=%.6f | CT=%.6f | CP=%.6f | eta=%.4f | T=%.3f N | P=%.3f W | AF=%.1f\n', ...
    J_isoP, CT_isoP, CP_isoP, eta_isoP, T_isoP, P_isoP, AF_apc);

fprintf('APC iso-tracao   : J=%.6f | CT=%.6f | CP=%.6f | eta=%.4f | T=%.3f N | P=%.3f W | AF=%.1f\n', ...
    J_isoT, CT_isoT, CP_isoT, eta_isoT, T_ref, P_isoT, AF_apc);

%% ===================== PLOTS ================================
figure('Color','w'); hold on; grid on;
plot(J_common, Cp4300, 'LineWidth', 1.5, 'DisplayName', 'APC C_P(J) @ 4300 rpm');
plot(J_same, CP_same, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Mesmo J');
plot(J_isoP, CP_isoP, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Iso-potencia');
plot(J_isoT, CP_isoT, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Iso-tracao');
xlabel('J'); ylabel('C_P');
title('APC 22x10 @ 4300 rpm');
legend('Location','best');

figure('Color','w'); hold on; grid on;
plot(J_common, T4300, 'LineWidth', 1.5, 'DisplayName', 'APC T(J) @ 4300 rpm');
plot(J_same, T_same, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Mesmo J');
plot(J_isoP, T_isoP, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Iso-potencia');
plot(J_isoT, T_ref,  'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Iso-tracao');
xlabel('J'); ylabel('T [N]');
title('APC 22x10 @ 4300 rpm');
legend('Location','best');

figure('Color','w'); hold on; grid on;
plot(J_common, P4300, 'LineWidth', 1.5, 'DisplayName', 'APC P(J) @ 4300 rpm');
plot(J_same, P_same, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Mesmo J');
plot(J_isoP, P_isoP, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Iso-potencia');
plot(J_isoT, P_isoT, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Iso-tracao');
xlabel('J'); ylabel('P [W]');
title('APC 22x10 @ 4300 rpm');
legend('Location','best');

figure('Color','w'); hold on; grid on;
plot(J_common, Pe4300, 'LineWidth', 1.5, 'DisplayName', 'APC \eta(J) @ 4300 rpm');
plot(J_same, eta_same, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Mesmo J');
plot(J_isoP, eta_isoP, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Iso-potencia');
plot(J_isoT, eta_isoT, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Iso-tracao');
xlabel('J'); ylabel('\eta');
title('APC 22x10 @ 4300 rpm');
legend('Location','best');
