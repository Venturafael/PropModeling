%==========================================================================
% MEL — projeto de hélice por Minimum Energy Loss (Adkins e Liebeck).
% A formulação, as constantes e a discretização do código do TCC foram
% mantidas. Apenas os caminhos locais foram organizados.
%
% Optional workspace inputs preserved across startup:
%   PE0_file, DAT_file, POL_file
%
% Original script description:
%   Adkins & Liebeck (MEL) — Otimização vs APC 22×10 (FIX + ξ)
%   - Polares Clark-Y (polares_clarky_expan.mat)
%   - Geometria APC (PE0) e dados experimentais (PER3)
%   - Resolve ζ para Pc alvo e varre Cl para AF >= AF_min
%   - Plota c(ξ), β(ξ), α(ξ), dCT/dξ, dCP/dξ e compara ponto final no plano J
%==========================================================================

clearvars -except PE0_file DAT_file POL_file
clc; close all

matlab_root = fileparts(mfilename('fullpath'));
addpath(fullfile(matlab_root, 'Functions'));

%% ===================== INPUTS (SI) =======================================
B        = 2;
D        = 0.559;         % [m] (22")
rho      = 1.225;         % [kg/m^3]
nu       = 1.47e-5;       % [m^2/s]
V        = 11.9;          % [m/s]
RPM      = 4300;          % [-]
n        = RPM/60;        % [1/s]

% alvo base
Cl_des_base = 0.700;
zeta_init   = 0.100;

% potência no eixo (alvo Pc) — Eq. L9b
P_axis = 0.50e3;          % [W]

% "afinadores"
alpha_bias_deg = -3.5;    % viés em α (β = φ + α)
Cd_factor      = 1.00;    % fator multiplicativo em Cd
c_min_tip      = 0.0;     % corda mínima na ponta [m]

% numérico / limites
maxit_zeta   = 60;  tol_zeta_abs = 1e-3;  tol_zeta_rel = 5e-4;
F_min        = 0.05;
a_limits     = [-0.20 0.60];
ap_limits    = [-0.50 0.50];
sphi_floor   = 1e-8;
eps_floor    = 1e-12;
a_sound      = 343;
tol_loc      = 1e-6;

% restrição de fabricabilidade (AF)
AF_min   = 70;
Cl_grid  = 0.50:0.02:1.40;

% arquivos
repo_paths = propeller_repo_paths();
if ~exist('PE0_file','var') || isempty(PE0_file)
    PE0_file = repo_paths.geometry_22x10;
end
if ~exist('DAT_file','var') || isempty(DAT_file)
    DAT_file = repo_paths.performance_22x10;
end
if ~exist('POL_file','var') || isempty(POL_file)
    POL_file = repo_paths.polar_clark_y;
end
assert(exist(PE0_file,'file')==2,'PE0 não encontrado');
assert(exist(DAT_file,'file')==2,'PER3 não encontrado');
assert(exist(POL_file,'file')==2,'Polares não encontrados');

%% ===================== POLARES CLARK-Y ==================================
Sload = load(POL_file);
if isfield(Sload,'ans') && isfield(Sload.ans,'perfil')
    pol = Sload.ans.perfil;
elseif isfield(Sload,'perfil')
    pol = Sload.perfil;
else
    error('Grid "perfil" não encontrado no MAT.');
end

%% ===================== GEOMETRIA APC (PE0) ===============================
APC = read_APC_PE0(PE0_file);
inch2m   = 0.0254;
R_apc_in = APC.RADIUS;                      % [in]
D_apc_m  = 2*R_apc_in*inch2m;               % [m]
xi_apc   = APC.STATION / R_apc_in;          % r/R
c_apc_m  = APC.CHORD   * inch2m;            % [m]
beta_apc_deg = atan( APC.PITCH_PRATHER ./ (2*pi*APC.STATION) ) * 180/pi;

AF_apc = (100000/16) * APC.BLADES * trapz(xi_apc, (c_apc_m/D_apc_m) .* (xi_apc.^3));
xi_apc_root = min(xi_apc);

%% ===================== MALHA RADIAL (alinhada) ===========================
xi_root = max(0.10, xi_apc_root);
xi      = [linspace(xi_root,0.95,14) 0.99]';
xi      = unique(min(max(xi, xi_root), 0.99));

%% ===================== DERIVADOS GLOBAIS ================================
R      = D/2;
J      = V/(n*D);
lambda = J/pi;
X      = xi./lambda;
Pc_target = 2*P_axis/(rho*V^3*pi*R^2);

P = struct('B',B,'D',D,'rho',rho,'nu',nu,'V',V,'n',n,'xi',xi,'R',R,'J',J, ...
           'lambda',lambda,'X',X,'alpha_bias_deg',alpha_bias_deg, ...
           'Cd_factor',Cd_factor,'F_min',F_min,'sphi_floor',sphi_floor, ...
           'eps_floor',eps_floor,'a_sound',a_sound,'a_limits',a_limits, ...
           'ap_limits',ap_limits,'tol_loc',tol_loc,'c_min_tip',c_min_tip);

%% ===================== ζ de base e varredura em Cl ======================
% ζ (base)
zeta_tmp = zeta_init;
for it = 1:maxit_zeta
    [~, CC] = estado_dado_zeta(zeta_tmp, Cl_des_base, pol, P);
    J1 = CC.J1;  J2 = max(CC.J2, eps);
    zeta_new = (-J1 + sqrt(J1^2 + 4*J2*Pc_target)) / (2*J2);
    if abs(zeta_new - zeta_tmp) < tol_zeta_abs || ...
       abs(zeta_new - zeta_tmp)/max(abs(zeta_new),1) < tol_zeta_rel
        zeta_tmp = zeta_new; break
    end
    zeta_tmp = zeta_new;
end
[S_base, ~] = estado_dado_zeta(zeta_tmp, Cl_des_base, pol, P);

%% ===================== VARREDURA EM Cl* — SALVANDO O HISTÓRICO ==========
% Cada Cl_i gera uma hélice completa: c(xi), beta(xi), alpha(xi), T, AF etc.
% Esse bloco guarda os escalares principais para mapear o espaço de projeto.

RES = repmat(struct( ...
    'Cl',[], 'T',[], 'ETA',[], 'AF',[], 'zeta',[], ...
    'CT',[], 'CP',[], 'Pc',[], ...
    'alpha_mean',[], 'alpha_min',[], 'alpha_max',[], ...
    'alpha_polar_mean',[], ...
    'c_mean',[], 'c_min',[], 'c_max',[]), numel(Cl_grid), 1);

for i = 1:numel(Cl_grid)

    Cl_i = Cl_grid(i);

    % Resolve zeta para cumprir a potência disponível
    zeta_loc = zeta_init;

    for it = 1:maxit_zeta

        [~, CC] = estado_dado_zeta(zeta_loc, Cl_i, pol, P);

        J1 = CC.J1;
        J2 = max(CC.J2, eps);

        zeta_new = (-J1 + sqrt(J1^2 + 4*J2*Pc_target)) / (2*J2);

        if abs(zeta_new - zeta_loc) < tol_zeta_abs || ...
           abs(zeta_new - zeta_loc)/max(abs(zeta_new),1) < tol_zeta_rel
            zeta_loc = zeta_new;
            break
        end

        zeta_loc = zeta_new;
    end

    % Estado completo da hélice para esse Cl_i
    [S_loc, ~] = estado_dado_zeta(zeta_loc, Cl_i, pol, P);

    % Cálculo do fator de atividade
    dxi          = gradient(S_loc.xi);
    AF_i_loc     = (100000/16) * B .* (S_loc.c./D) .* (S_loc.xi.^3) .* dxi;
    AF_tot_loc   = sum(AF_i_loc,'omitnan');

    % Estatísticas de alpha
    % ATENÇÃO:
    % Se S_loc.alpha já inclui alpha_bias_deg, então alpha_polar = alpha - alpha_bias_deg.
    % Como alpha_bias_deg = -3.5, isso desloca o alpha "real da polar" para cima.
    if isfield(S_loc,'alpha')
        alpha_vec = S_loc.alpha(:);

        alpha_mean_loc = mean(alpha_vec,'omitnan');
        alpha_min_loc  = min(alpha_vec,[],'omitnan');
        alpha_max_loc  = max(alpha_vec,[],'omitnan');

        alpha_polar_mean_loc = mean(alpha_vec - alpha_bias_deg,'omitnan');
    else
        alpha_mean_loc = NaN;
        alpha_min_loc  = NaN;
        alpha_max_loc  = NaN;
        alpha_polar_mean_loc = NaN;
    end

    % Estatísticas de corda
    c_mean_loc = mean(S_loc.c,'omitnan');
    c_min_loc  = min(S_loc.c,[],'omitnan');
    c_max_loc  = max(S_loc.c,[],'omitnan');

    % Salva resultados escalares
    RES(i).Cl   = Cl_i;
    RES(i).T    = S_loc.T;
    RES(i).ETA  = S_loc.ETA;
    RES(i).AF   = AF_tot_loc;
    RES(i).zeta = S_loc.zeta_tip;
    RES(i).CT   = S_loc.CT;
    RES(i).CP   = S_loc.CP;
    RES(i).Pc   = S_loc.PC;

    RES(i).alpha_mean       = alpha_mean_loc;
    RES(i).alpha_min        = alpha_min_loc;
    RES(i).alpha_max        = alpha_max_loc;
    RES(i).alpha_polar_mean = alpha_polar_mean_loc;

    RES(i).c_mean = c_mean_loc;
    RES(i).c_min  = c_min_loc;
    RES(i).c_max  = c_max_loc;
end

feas = [RES.AF] >= AF_min;
if any(feas)
    feas_idx   = find(feas);
    [~,kT]     = max([RES(feas).T]);
    idx_bestT  = feas_idx(kT);
else
    [~, idx_bestT] = max([RES.AF]);   % fallback
    warning('Nenhum caso com AF >= %.1f. Reportando o de maior AF.', AF_min);
end
Cl_star = RES(idx_bestT).Cl;
if ~isscalar(Cl_star) || ~isfinite(Cl_star), Cl_star = Cl_des_base; end

% resolve FINAL com Cl*
zeta = zeta_init;
for it = 1:maxit_zeta
    [~, CC] = estado_dado_zeta(zeta, Cl_star, pol, P);
    J1 = CC.J1;  J2 = max(CC.J2, eps);
    zeta_new = (-J1 + sqrt(J1^2 + 4*J2*Pc_target)) / (2*J2);
    if abs(zeta_new - zeta) < tol_zeta_abs || ...
       abs(zeta_new - zeta)/max(abs(zeta_new),1) < tol_zeta_rel
        zeta = zeta_new; break
    end
    zeta = zeta_new;
end
[S_final, ~] = estado_dado_zeta(zeta, Cl_star, pol, P);

% AF FINAL
dxi_final      = gradient(S_final.xi);
AF_i_final     = (100000/16) * B .* (S_final.c./D) .* (S_final.xi.^3) .* dxi_final;
AF_total_final = sum(AF_i_final,'omitnan');

%% ===================== Dados experimentais (PER3, RPM mais próximo) =====
EXP = read_APC_PER3_nearest(DAT_file, RPM);
J_exp = EXP.J; Pe_exp = EXP.Pe; Ct_exp = EXP.Ct; Cp_exp = EXP.Cp; target_rpm = EXP.target_rpm;

% nosso ponto
J_final  = S_final.J;
Ct_final = S_final.CT;
Cp_final = S_final.CP;
Pe_final = S_final.ETA;

% APC no mesmo J (interp)
Jq        = min(max(J_final, min(J_exp)), max(J_exp));
Ct_apc_q  = interp1(J_exp, Ct_exp, Jq, 'pchip');
Cp_apc_q  = interp1(J_exp, Cp_exp, Jq, 'pchip');
Pe_apc_q  = interp1(J_exp, Pe_exp, Jq, 'pchip');

n_apc  = target_rpm/60;
T_apc_q = Ct_apc_q * rho * n_apc^2 * D_apc_m^4;
T_fin   = S_final.T;

%% ===================== Estado MEL com GEOMETRIA da APC ===================
P_apc = P; P_apc.xi = xi_apc(:); P_apc.X = P_apc.xi ./ P_apc.lambda;

zeta_apc = zeta_init;
for it = 1:maxit_zeta
    [~, CCa] = estado_dado_zeta_geom(zeta_apc, pol, P_apc, beta_apc_deg(:), c_apc_m(:));
    J1a = CCa.J1;  J2a = max(CCa.J2, eps);
    zeta_new = (-J1a + sqrt(J1a^2 + 4*J2a*Pc_target)) / (2*J2a);
    if abs(zeta_new - zeta_apc) < tol_zeta_abs || ...
       abs(zeta_new - zeta_apc)/max(abs(zeta_new),1) < tol_zeta_rel
        zeta_apc = zeta_new; break
    end
    zeta_apc = zeta_new;
end
[S_apc, ~] = estado_dado_zeta_geom(zeta_apc, pol, P_apc, beta_apc_deg(:), c_apc_m(:));

%% ===================== PRINT RESUMO =====================================
fprintf('\n===== RESUMO (FINAL vs APC) =====\n');
fprintf('Nosso FINAL:   Cl*=%.3f | ζ=%.4f | J=%.3f | CT=%.4f | CP=%.4f | η=%.3f | T=%.2f N | AF=%.1f\n', ...
        Cl_star, S_final.zeta_tip, J_final, Ct_final, Cp_final, Pe_final, T_fin, AF_total_final);
fprintf('APC 22x10:     B=%g | D=%.3f m | AF≈%.1f | início de perfil ≈ %.3fR\n', APC.BLADES, D_apc_m, AF_apc, xi_apc_root);
fprintf('APC (interp.): J=%.3f | CT=%.4f | CP=%.4f | η=%.3f | T(%.0f RPM)=%.2f N\n', ...
        Jq, Ct_apc_q, Cp_apc_q, Pe_apc_q, target_rpm, T_apc_q);

%% ===================== PLOT STYLE =======================================
fontName = 'Helvetica'; fs_axis = 12; fs_title = 14;
axlw = 1.0; gridA = 0.18; minorA = 0.12; figSize_cm = [14 9];
newfig = @(name) set(figure('Name',name,'Color','w','Units','centimeters', ...
                     'Position',[1 1 figSize_cm]),'PaperPositionMode','auto');
format_axes = @(ax) set(ax,'FontName',fontName,'FontSize',fs_axis,'Box','off', ...
                     'TickDir','out','LineWidth',axlw,'GridAlpha',gridA,'MinorGridAlpha',minorA, ...
                     'XMinorTick','on','YMinorTick','on');
c_final = [0.20 0.60 0.40];
c_apc   = [0.10 0.25 0.55];

%% ===================== MAPAS DA VARREDURA Cl* × AF × T ===================

% Vetores da varredura
Cl_vec    = [RES.Cl];
T_vec     = [RES.T];
AF_vec    = [RES.AF];
ETA_vec   = [RES.ETA];
CT_vec    = [RES.CT];
CP_vec    = [RES.CP];
zeta_vec  = [RES.zeta];

alpha_mean_vec       = [RES.alpha_mean];
alpha_min_vec        = [RES.alpha_min];
alpha_max_vec        = [RES.alpha_max];
alpha_polar_mean_vec = [RES.alpha_polar_mean];

feas_vec = AF_vec >= AF_min;

T_star   = RES(idx_bestT).T;
AF_star  = RES(idx_bestT).AF;
ETA_star = RES(idx_bestT).ETA;


%% ------------------------------------------------------------------------
%  GRÁFICO 1 — T e AF em função de Cl*
%  Esse é o gráfico mais importante para entender por que o código escolheu Cl*=0.76
% -------------------------------------------------------------------------

newfig('Mapa 2D — Tração e AF por Cl*'); hold on; grid on

yyaxis left
plot(Cl_vec, T_vec, '-o', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Tração T');

plot(Cl_vec(feas_vec), T_vec(feas_vec), 'o', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Casos viáveis: AF \geq AF_{min}');

plot(Cl_star, T_star, 'p', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf('Ótimo: C_l^*=%.3f', Cl_star));

ylabel('Tração T [N]');

yyaxis right
plot(Cl_vec, AF_vec, '-s', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Fator de Atividade AF');

yline(AF_min, '--', ...
    sprintf('AF_{min}=%.0f', AF_min), ...
    'LineWidth', 1.3, ...
    'DisplayName', 'AF mínimo');

ylabel('Fator de atividade AF [-]');

xline(Cl_star, ':', ...
    sprintf('C_l^*=%.3f', Cl_star), ...
    'LineWidth', 1.4, ...
    'DisplayName', 'C_l^* escolhido');

xlabel('C_l^* desejado [-]');
title('Varredura do projeto MEL: tração e AF em função de C_l^*', ...
      'FontSize', fs_title, 'FontName', fontName);

legend('Location','best','Box','off');
format_axes(gca);


%% ------------------------------------------------------------------------
%  GRÁFICO 2 — Espaço 3D: Cl*, AF e T
%  Observação: como AF é consequência do Cl*, isso é uma curva 3D, não uma superfície física.
% -------------------------------------------------------------------------

newfig('Espaço 3D — Cl* × AF × T'); hold on; grid on

% Curva paramétrica no espaço de projeto
plot3(Cl_vec, AF_vec, T_vec, '-', ...
    'LineWidth', 1.2, ...
    'Color', [0.25 0.25 0.25], ...
    'DisplayName', 'Caminho da varredura');

% Nuvem colorida pela eficiência
scatter3(Cl_vec, AF_vec, T_vec, 55, ETA_vec, ...
    'filled', ...
    'DisplayName', 'Casos avaliados');

% Pontos viáveis destacados
scatter3(Cl_vec(feas_vec), AF_vec(feas_vec), T_vec(feas_vec), 80, ETA_vec(feas_vec), ...
    'o', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Viáveis: AF \geq AF_{min}');

% Ótimo escolhido
plot3(Cl_star, AF_star, T_star, 'p', ...
    'MarkerSize', 16, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf('Ótimo: C_l^*=%.3f', Cl_star));

% Plano visual da restrição AF_min
[CL_plane, T_plane] = meshgrid( ...
    linspace(min(Cl_vec), max(Cl_vec), 2), ...
    linspace(min(T_vec), max(T_vec), 2));

AF_plane = AF_min * ones(size(CL_plane));

surf(CL_plane, AF_plane, T_plane, ...
    'FaceAlpha', 0.10, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Plano AF_{min}');

xlabel('C_l^* desejado [-]');
ylabel('Fator de atividade AF [-]');
zlabel('Tração T [N]');
title('Espaço de projeto: C_l^*, AF e tração', ...
      'FontSize', fs_title, 'FontName', fontName);

cb = colorbar;
ylabel(cb, '\eta [-]');

view(45, 25);
legend('Location','best','Box','off');
format_axes(gca);


%% ------------------------------------------------------------------------
%  GRÁFICO 3 — Bônus: alpha médio/min/max em função de Cl*
%  Esse gráfico ajuda a responder por que o alpha ficou tão baixo.
% -------------------------------------------------------------------------

if any(isfinite(alpha_mean_vec))

    newfig('Alpha operacional por Cl*'); hold on; grid on

    % Faixa radial de alpha: min até max ao longo da pá
    fill([Cl_vec fliplr(Cl_vec)], ...
         [alpha_min_vec fliplr(alpha_max_vec)], ...
         [0.80 0.80 0.80], ...
         'FaceAlpha', 0.35, ...
         'EdgeColor', 'none', ...
         'DisplayName', '\alpha_{min} até \alpha_{max}');

    plot(Cl_vec, alpha_mean_vec, '-o', ...
        'LineWidth', 1.8, ...
        'DisplayName', '\alpha médio usado na geometria');

    plot(Cl_vec, alpha_polar_mean_vec, '--s', ...
        'LineWidth', 1.6, ...
        'DisplayName', '\alpha médio estimado da polar');

    plot(Cl_star, RES(idx_bestT).alpha_mean, 'p', ...
        'MarkerSize', 14, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', sprintf('Ótimo: C_l^*=%.3f', Cl_star));

    yline(0, ':', '\alpha = 0^\circ', ...
        'LineWidth', 1.2, ...
        'DisplayName', '\alpha = 0^\circ');

    xline(Cl_star, ':', ...
        sprintf('C_l^*=%.3f', Cl_star), ...
        'LineWidth', 1.4, ...
        'DisplayName', 'C_l^* escolhido');

    xlabel('C_l^* desejado [-]');
    ylabel('\alpha [deg]');
    title('\alpha operacional ao longo da varredura de C_l^*', ...
          'FontSize', fs_title, 'FontName', fontName);

    legend('Location','best','Box','off');
    format_axes(gca);
end


%% ------------------------------------------------------------------------
%  GRÁFICO 4 — Bônus: CT, CP e eta em função de Cl*
%  Ajuda a enxergar se a tração cresce, mas a potência/arrasto penaliza.
% -------------------------------------------------------------------------

newfig('CT CP eta por Cl*');
t2 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on
plot(Cl_vec, CT_vec, '-o', 'LineWidth', 1.6, 'DisplayName', 'C_T');
xline(Cl_star, ':', sprintf('C_l^*=%.3f', Cl_star), 'LineWidth', 1.2);
ylabel('C_T [-]');
legend('Location','best','Box','off');
format_axes(gca);

nexttile; hold on; grid on
plot(Cl_vec, CP_vec, '-s', 'LineWidth', 1.6, 'DisplayName', 'C_P');
xline(Cl_star, ':', sprintf('C_l^*=%.3f', Cl_star), 'LineWidth', 1.2);
ylabel('C_P [-]');
legend('Location','best','Box','off');
format_axes(gca);

nexttile; hold on; grid on
plot(Cl_vec, ETA_vec, '-^', 'LineWidth', 1.6, 'DisplayName', '\eta');
xline(Cl_star, ':', sprintf('C_l^*=%.3f', Cl_star), 'LineWidth', 1.2);
xlabel('C_l^* desejado [-]');
ylabel('\eta [-]');
legend('Location','best','Box','off');
format_axes(gca);

title(t2, 'Coeficientes globais ao longo da varredura em C_l^*', ...
      'FontSize', fs_title, 'FontName', fontName);

%% ===================== GRÁFICOS — GEOMETRIA ==============================
% c(ξ)
newfig('c(\xi): FINAL vs APC'); hold on; grid on
plot(S_final.xi, S_final.c, 'o-','LineWidth',1.8,'Color',c_final,'DisplayName','FINAL (ótimo)');
plot(xi_apc,     c_apc_m,   's-','LineWidth',1.6,'Color',c_apc  ,'DisplayName','APC 22×10');
xlabel('\xi = r/R'); ylabel('c [m]');
title(sprintf('Distribuição de corda — AF_{final}=%.1f | AF_{APC}=%.1f', AF_total_final, AF_apc), ...
      'FontSize',fs_title,'FontName',fontName);
legend('Location','best','Box','off'); format_axes(gca);

% β(ξ)
newfig('\beta(\xi): FINAL vs APC'); hold on; grid on
plot(S_final.xi, S_final.beta, 'o-','LineWidth',1.8,'Color',c_final,'DisplayName','FINAL (\beta=\phi+\alpha)');
plot(xi_apc,     beta_apc_deg,'s-','LineWidth',1.6,'Color',c_apc  ,'DisplayName','APC 22×10 (via pitch)');
xlabel('\xi = r/R'); ylabel('\beta [deg]');
title('\beta(\xi) — comparação geométrica','FontSize',fs_title,'FontName',fontName);
legend('Location','best','Box','off'); format_axes(gca);

% α(ξ)
newfig('\alpha(\xi): FINAL'); hold on; grid on
plot(S_final.xi, S_final.alpha, 'd-','LineWidth',1.8,'Color',c_final,'DisplayName','FINAL');
xlabel('\xi = r/R'); ylabel('\alpha [deg]');
title('\alpha(\xi) — caso FINAL','FontSize',fs_title,'FontName',fontName);
legend('Location','best','Box','off'); format_axes(gca);

%% ===================== GRÁFICOS — DISTRIBUIÇÕES =========================
% dCT/dξ (nosso)
newfig('dCT/d\xi — nossa hélice'); hold on; grid on
plot(S_final.xi, S_final.dCTdxi, '-o','LineWidth',1.8,'Color',c_final);
xlabel('\xi = r/R'); ylabel('dC_T/d\xi [-]'); title('Distribuição radial de C_T (nosso)','FontSize',fs_title);
format_axes(gca);

% dCP/dξ (nosso)
newfig('dCP/d\xi — nossa hélice'); hold on; grid on
plot(S_final.xi, S_final.dCPdxi, '-s','LineWidth',1.8,'Color',c_final);
xlabel('\xi = r/R'); ylabel('dC_P/d\xi [-]'); title('Distribuição radial de C_P (nosso)','FontSize',fs_title);
format_axes(gca);

% dCT/dξ (comparação)
newfig('dCT/d\xi — nosso vs APC'); hold on; grid on
plot(S_final.xi, S_final.dCTdxi, '-o','LineWidth',1.8,'Color',c_final,'DisplayName','Nosso');
plot(S_apc.xi ,  S_apc.dCTdxi ,  '-s','LineWidth',1.4,'Color',c_apc  ,'DisplayName','APC (MEL+geom)');
xlabel('\xi = r/R'); ylabel('dC_T/d\xi [-]'); title('Distribuição radial de C_T — comparação','FontSize',fs_title);
legend('Location','best','Box','off'); format_axes(gca);

% dCP/dξ (comparação)
newfig('dCP/d\xi — nosso vs APC'); hold on; grid on
plot(S_final.xi, S_final.dCPdxi, '-o','LineWidth',1.8,'Color',c_final,'DisplayName','Nosso');
plot(S_apc.xi ,  S_apc.dCPdxi ,  '-s','LineWidth',1.4,'Color',c_apc  ,'DisplayName','APC (MEL+geom)');
xlabel('\xi = r/R'); ylabel('dC_P/d\xi [-]'); title('Distribuição radial de C_P — comparação','FontSize',fs_title);
legend('Location','best','Box','off'); format_axes(gca);

%% ===================== GRÁFICOS — EXPERIMENTAL (vs J) ===================
newfig('Curvas experimentais APC + ponto FINAL');
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on
plot(J_exp, Ct_exp,'-o','LineWidth',1.2,'Color',[0.1 0.25 0.55],'DisplayName','APC C_T');
plot(J_final, Ct_final,'p','MarkerSize',10,'MarkerFaceColor',[0.7 0.15 0.15], ...
     'MarkerEdgeColor',[0.7 0.15 0.15],'DisplayName','FINAL (nosso)');
ylabel('C_T'); legend('Location','best','Box','off'); format_axes(gca);

nexttile; hold on; grid on
plot(J_exp, Cp_exp,'-s','LineWidth',1.2,'Color',[0.40 0.00 0.40],'DisplayName','APC C_P');
plot(J_final, Cp_final,'p','MarkerSize',10,'MarkerFaceColor',[0.7 0.15 0.15], ...
     'MarkerEdgeColor',[0.7 0.15 0.15],'DisplayName','FINAL (nosso)');
ylabel('C_P'); legend('Location','best','Box','off'); format_axes(gca);

nexttile; hold on; grid on
plot(J_exp, Pe_exp,'-^','LineWidth',1.2,'Color',[0 0.6 0.2],'DisplayName','APC \eta');
plot(J_final, Pe_final,'p','MarkerSize',10,'MarkerFaceColor',[0.7 0.15 0.15], ...
     'MarkerEdgeColor',[0.7 0.15 0.15],'DisplayName','FINAL (nosso)');
xlabel('J'); ylabel('\eta'); legend('Location','best','Box','off'); format_axes(gca);

title(t, sprintf('APC 22×10 @ %d RPM (interp) × Ponto Ótimo (J=%.3f)', target_rpm, J_final), ...
      'FontSize',fs_title,'FontName',fontName);
