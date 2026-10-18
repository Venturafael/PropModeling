%======================================================================
% VORTEX — modelo de McCormick usado no TCC.
% A formulação, a discretização e as constantes do código original foram
% mantidas. Apenas os caminhos locais foram organizados.
%
% O teste histórico de espessura foi mantido. Com o leitor fornecido, o
% código usa o fallback t/c = 0.103, como na versão do TCC.
%
% Expected optional workspace inputs:
%   geometry_file, performance_file, polar_file, target_rpm, output_dir
%
% Original script description:
% (McCormick + Vortex; usa polares Clark-Y via busca_polares_corretas)
%   1) Lê geometria e dados de referência
%   2) Calcula mapa teórico com McCormick/Lifting-Line [núcleo inalterado]
%   3) Métricas globais (R², r² e sMAPE)  [APENAS NA TABELA]
%   4) Mapas (Teórico, Referência, Teo+ηref) e comparativos (sem sMAPE)
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

%% 1) GEOMETRIA -----------------------------------------------------------
run(fullfile(repo_paths.functions, 'leitura_geometria_pe0.m')); % STATION [in], CHORD [in], PITCH_LE_TE [in], RADIUS [in], BLADES

B = BLADES;
R_in = RADIUS;                 % [in]
D_in = 2*R_in;                 % [in]
deg = pi/180;

% garantir vetores-coluna
STATION = STATION(:);
CHORD = CHORD(:);
PITCH_LE_TE = PITCH_LE_TE(:);

x_rel = STATION./R_in;         % r/R
c_in = CHORD;                  % [in]
pitch = PITCH_LE_TE;           % [in]
n_x = numel(x_rel);

% espessura t/c se disponível; fallback 0.103
if     exist('T_OVER_C','var'), tc_vec = T_OVER_C(:);
elseif exist('TC','var'),       tc_vec = TC(:);
elseif exist('t_over_c','var'), tc_vec = t_over_c(:);
else,                           tc_vec = 0.103*ones(n_x,1);
end

% geometria em ângulo
beta_deg = atan2(pitch, 2*pi*STATION)/deg;

% tip data para Prandtl (φ fixo na ponta)
[~, i_tip] = max(x_rel);
i_tc_tip = max(i_tip-1, 1);                  % usa ~0.99R se existir
beta_tip_deg = beta_deg(i_tip);
alpha0_tip = 41*tc_vec(i_tc_tip);
phi_t_corr = (beta_tip_deg + alpha0_tip)*deg; % φ_t,corr (apenas p/ Prandtl)
sin_phi_tcorr = sin(phi_t_corr);

%% 2) DADOS DE REFERÊNCIA (FABRICANTE/EXPER.) ----------------------------
run(fullfile(repo_paths.functions, 'leitura_desempenho_dat.m')); % J_exp, Ct_exp, Cp_exp, Pe_exp, target_rpm

J_ref = J_exp(:)';              % mantém os mesmos J do ensaio
Ct_ref = Ct_exp(:);
Cp_ref = Cp_exp(:);
Eta_ref = Pe_exp(:);
n_J = numel(J_ref);

%% 3) POLARES + PROPRIEDADES DO AR ---------------------------------------
% Usa o banco no formato ans.perfil (NR x NM)
S = load(polar_file);
pol = S.ans.perfil;

rho = 1.225;   % kg/m^3
mu = 1.81e-5;  % Pa·s
Vson = 343;    % m/s

%% 4) ARRAYS DE SAÍDA + DIAGNÓSTICOS -------------------------------------
CT_the = zeros(n_J,1);
CP_the = zeros(n_J,1);
ETA_the = NaN(n_J,1);           % só preenche quando CT>0 & CP>0

CL_x = NaN(n_x, n_J);
CD_x = NaN(n_x, n_J);
ALF_x = NaN(n_x, n_J);          % α "usado" na polar
RE_x = NaN(n_x, n_J);
MA_x = NaN(n_x, n_J);

%% 5) LOOP EM J  [núcleo McCormick/Lifting-Line INALTERADO] --------------
for j = 1:n_J
    J = J_ref(j);
    lam = J/pi;

    dCTdx = zeros(n_x,1);
    dCPdx = zeros(n_x,1);

    for k = 1:n_x
        x = x_rel(k);
        c = c_in(k);
        tc = tc_vec(k);
        bet = beta_deg(k);

        % razão de área média (modelo McCormick original)
        sigma_m = B*c/(pi*R_in);

        % iteração ponto-fixo para wt/VT
        wt = 0.010;
        convergiu = false;
        VE_VT = NaN; wa = NaN; phi_rel = NaN; phi_geo = NaN;

        for it = 1:500
            % V1: wa/VT
            radic = lam^2 + 4*wt*(x - wt);
            if radic <= 0, break; end
            wa = 0.5*(-lam + sqrt(radic));

            % ângulos (em rad)
            phi_geo = atan(lam/x);
            phi_rel = atan2((lam + wa), (x - wt));

            % correções aerodinâmicas (em graus)
            Delta_theta = (atan((lam+wa)/(x-2*wt)) - atan((lam+wa)/x))/deg;
            Delta_alpha = (4*lam*sigma_m/(15*(lam^2 + x^2)))*tc*(180/pi);

            % (1) Ângulo para a polar — geométrico
            alpha_geo = bet - rad2deg(phi_rel) - Delta_alpha;

            % Re/Mach locais
            [Re_i_raw, Ma_i_raw] = reynolds_mach_local( ...
                STATION(k), c, J, target_rpm, R_in, rho, mu, Vson);

            % (2) Lookup das polares (robusto + monotonia)
            [Cl_p, Cd, ~, a_use, Re_u, Ma_u] = busca_polares_corretas( ...
                alpha_geo, Re_i_raw, Ma_i_raw, pol, 'linear', k, j, true, true);

            % (3) Correção McCormick por Δθ na própria polar
            a_loc = dcl_dalpha_local(alpha_geo, Re_u, Ma_u, pol, 'linear', k, j, true);
            if ~isfinite(a_loc) || a_loc <= 0
                a_loc = 0.1*(1+tc);   % fallback suave
            end
            Cl = Cl_p - a_loc*(Delta_theta/4);

            % módulo de velocidade relativa adim.
            VE_VT = hypot(lam+wa, x-wt);

            % Prandtl tip-φ (constante em x; usa φ_t,corr do tip)
            f = (B/2)*(1 - x)/max(sin_phi_tcorr, 1e-8);
            F = max(2/pi*acos(exp(-max(f,0))), 0.05);

            % ponto-fixo wt novo (limitador simétrico + sub-relaxação)
            wt_new = sigma_m*Cl*VE_VT/(8*x*F);
            wt_new = max(min(wt_new, 0.9*x), -0.9*x);
            wt_new = 0.55*wt_new + 0.45*wt;

            if abs(wt_new - wt) < 1e-5
                wt = wt_new; convergiu = true; break;
            end
            wt = wt_new;

            % diagnósticos do passo
            ALF_x(k,j) = a_use;
            RE_x(k,j) = Re_u;
            MA_x(k,j) = Ma_u;
            CL_x(k,j) = Cl;
            CD_x(k,j) = Cd;
        end

        if ~convergiu || ~isfinite(Cl) || ~isfinite(Cd) || ~isfinite(phi_rel)
            dCTdx(k) = 0; dCPdx(k) = 0;
            if ~isfinite(ALF_x(k,j)), ALF_x(k,j) = alpha_geo; end
            if ~isfinite(RE_x(k,j)),  RE_x(k,j)  = Re_i_raw;  end
            if ~isfinite(MA_x(k,j)),  MA_x(k,j)  = Ma_i_raw;  end
            continue
        end

        % salva diagnósticos finais do ponto (se ainda faltou)
        if ~isfinite(ALF_x(k,j)), ALF_x(k,j) = a_use; end
        if ~isfinite(RE_x(k,j)),  RE_x(k,j)  = Re_u;  end
        if ~isfinite(MA_x(k,j)),  MA_x(k,j)  = Ma_u;  end
        if ~isfinite(CL_x(k,j)),  CL_x(k,j)  = Cl;    end
        if ~isfinite(CD_x(k,j)),  CD_x(k,j)  = Cd;    end

        % incrementos (V22–V23)
        J2pi2 = J^2 + (pi*x)^2;
        dCTdx(k) = (pi/8)*J2pi2*sigma_m*(Cl*cos(phi_rel) - Cd*sin(phi_rel));
        dCPdx(k) = (pi/8)*(pi*x)*J2pi2*sigma_m*(Cl*sin(phi_rel) + Cd*cos(phi_rel));
    end

    % interpola NaNs ao longo de x antes de integrar
    dCTdx = linfill1d(dCTdx);
    dCPdx = linfill1d(dCPdx);
    ALF_x(:,j) = linfill1d(ALF_x(:,j));
    RE_x(:,j) = linfill1d(RE_x(:,j));
    MA_x(:,j) = linfill1d(MA_x(:,j));
    CL_x(:,j) = linfill1d(CL_x(:,j));
    CD_x(:,j) = linfill1d(CD_x(:,j));

    % integração radial
    CT_the(j) = trapz(x_rel, dCTdx);
    CP_the(j) = trapz(x_rel, dCPdx);

    % eficiência somente em regime propulsivo
    if (CT_the(j) > 0) && (CP_the(j) > 0)
        ETA_the(j) = J*CT_the(j)/CP_the(j);
    end
end

% interpola NaNs ao longo de J (diagnósticos Re/Mach)
for i = 1:n_x
    RE_x(i,:)  = linfill1d(RE_x(i,:).').';
    MA_x(i,:)  = linfill1d(MA_x(i,:).').';
    ALF_x(i,:) = linfill1d(ALF_x(i,:).').';
    CL_x(i,:)  = linfill1d(CL_x(i,:).').';
    CD_x(i,:)  = linfill1d(CD_x(i,:).').';
end

%% 6) MÉTRICAS (R², r² e sMAPE) — APENAS NA TABELA -----------------------
ETA_mask = ETA_the;
bad = ~(CT_the > 0 & CP_the > 0);
ETA_mask(bad) = NaN;

met_ct = series_metrics_base(Ct_ref, CT_the);
met_cp = series_metrics_base(Cp_ref, CP_the);
met_eta = series_metrics_base(Eta_ref, ETA_mask);

%% 7) PLOTS — parâmetros visuais (padrão Compare_*) ----------------------
fontName = 'Helvetica'; fs_axis = 12; fs_title = 14;
axlw = 1.0; gridA = 0.18; minorA = 0.12; figSize_cm = [14 9];

lw_ct = 2.8; lw_cp = 2.8; lw_eta = 2.8; lw_ref = 2.8;
ms_ref = 3.2; mk_ref = 'o';

ls_ct = '-'; ls_cp = '--'; ls_eta = '-.'; ls_ref = ':';

co_map = [0.20 0.60 0.40; 0.70 0.15 0.15; 0.10 0.25 0.55];
cp_color = [0.40 0.00 0.40];

cmap = @(k) co_map(mod(k-1, size(co_map,1))+1, :);
lighten = @(rgb,a) (1-a)*rgb + a*[1 1 1];

format_axes = @(ax) set(ax,'FontName',fontName,'FontSize',fs_axis,'Box','off', ...
    'TickDir','out','LineWidth',axlw,'GridAlpha',gridA,'MinorGridAlpha',minorA, ...
    'XMinorTick','on','YMinorTick','on');
newfig = @(name) set(figure('Name',name,'Color','w','Units','centimeters', ...
    'Position',[1 1 figSize_cm]),'PaperPositionMode','auto');
padlims = @(v,lo,hi) local_padlims(v,lo,hi);
leg = {'Location','best','Box','off'};

lab_CT_ref = 'C_T (ref)';
lab_CP_ref = 'C_P (ref)';
lab_ETA_ref = '\eta (ref)';

% ------------------- FIGURAS DE MAPA ------------------------------------
% [Figura 1] Mapa Teórico
newfig('Mapa Teórico'); hold on; grid on
yyaxis left
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cmap(2), 'DisplayName', 'C_P (teo)');
ylabel('C_T , C_P'); ylim(padlims([CT_the; CP_the; 0], 0.05, 0.08))
yyaxis right
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
ylabel('Eficiência \eta'); ylim([0 1])
xlabel('Avanço J')
title(sprintf('Mapa Teórico – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName)
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = cmap(3); ax.YAxis(1).Color = [0.2 0.2 0.2];

% [Figura 2] Mapa de Referência (fabricante/experimental)
newfig('Mapa de Referência'); hold on; grid on
yyaxis left
plot(J_ref, Ct_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', cmap(1), 'MarkerEdgeColor', cmap(1), 'Color', lighten(cmap(1), 0.2), 'DisplayName', lab_CT_ref);
plot(J_ref, Cp_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', cp_color, 'MarkerEdgeColor', cp_color, 'Color', lighten(cp_color, 0.2), 'DisplayName', lab_CP_ref);
ylabel('C_T , C_P'); ylim(padlims([Ct_ref; Cp_ref; 0], 0.05, 0.08))
yyaxis right
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('Eficiência \eta'); ylim([0 1])
xlabel('Avanço J')
title(sprintf('Mapa de Referência – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName)
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = [0 0 0]; ax.YAxis(1).Color = [0.2 0.2 0.2];

% [Figura 3] Mapa Teórico + η (ref)
newfig('Mapa Teo + \eta (ref)'); hold on; grid on
yyaxis left
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cmap(2), 'DisplayName', 'C_P (teo)');
ylabel('C_T , C_P'); ylim(padlims([CT_the; CP_the; 0], 0.05, 0.08))
yyaxis right
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('Eficiência \eta'); ylim([0 1])
xlabel('Avanço J')
title(sprintf('Mapa Teórico + \\eta (ref) – RPM %d', target_rpm), 'FontSize', fs_title, 'FontName', fontName)
legend(leg{:}); ax = gca; format_axes(ax); ax.YAxis(2).Color = [0 0 0]; ax.YAxis(1).Color = [0.2 0.2 0.2];

% ------------------- COMPARATIVOS (sem sMAPE/eixo direito) --------------
% [Figura 4] C_T × J
newfig('C_T × J – Teo × Ref'); hold on; grid on
plot(J_ref, CT_the, ls_ct, 'LineWidth', lw_ct, 'Color', cmap(1), 'DisplayName', 'C_T (teo)');
plot(J_ref, Ct_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', cmap(1), 'MarkerEdgeColor', cmap(1), 'Color', lighten(cmap(1), 0.2), 'DisplayName', lab_CT_ref);
ylabel('C_T'); ylim(padlims([CT_the; Ct_ref; 0], 0.05, 0.08))
xlabel('Avanço J'); title(sprintf('C_T × J – Teo × Ref  (RPM %d)', target_rpm))
legend(leg{:}); format_axes(gca)

% [Figura 5] C_P × J
newfig('C_P × J – Teo × Ref'); hold on; grid on
plot(J_ref, CP_the, ls_cp, 'LineWidth', lw_cp, 'Color', cp_color, 'DisplayName', 'C_P (teo)');
plot(J_ref, Cp_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', cp_color, 'MarkerEdgeColor', cp_color, 'Color', lighten(cp_color, 0.2), 'DisplayName', lab_CP_ref);
ylabel('C_P'); ylim(padlims([CP_the; Cp_ref; 0], 0.05, 0.08))
xlabel('Avanço J'); title(sprintf('C_P × J – Teo × Ref  (RPM %d)', target_rpm))
legend(leg{:}); format_axes(gca)

% [Figura 6] η × J
newfig('\eta × J – Teo × Ref'); hold on; grid on
plot(J_ref, ETA_the, ls_eta, 'LineWidth', lw_eta, 'Color', cmap(3), 'DisplayName', '\eta (teo)');
plot(J_ref, Eta_ref, ls_ref, 'LineWidth', lw_ref, 'Marker', mk_ref, 'MarkerSize', ms_ref, ...
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'DisplayName', lab_ETA_ref);
ylabel('\eta'); ylim([0 1])
xlabel('Avanço J'); title(sprintf('\\eta × J – Teo × Ref  (RPM %d)', target_rpm))
legend(leg{:}); format_axes(gca)

%% 8) TABELA RESUMO (console + CSV opcional) ------------------------------
labels = {'C_T','C_P','\eta'};
csvname = fullfile(output_dir, sprintf('metrics_vortex_RPM_%d.csv', target_rpm));

smape_ct_final = met_ct.sMAPE_acc_final;
smape_cp_final = met_cp.sMAPE_acc_final;
smape_eta_final = met_eta.sMAPE_acc_final;

print_metrics_table_acc( ...
    met_ct, met_cp, met_eta, ...
    [smape_ct_final smape_cp_final smape_eta_final], ...
    labels, target_rpm, csvname);

%% 9) DIAGNÓSTICOS RADIAIS (Re/Mach) -------------------------------------
co_rad = [0.20 0.60 0.40; 0.20 0.40 0.70; 0.50 0.40 0.70; 0.65 0.45 0.25; 0.90 0.70 0.30];
crad = @(k) co_rad(mod(k-1,size(co_rad,1))+1,:);

nJshow = min(6, n_J);  nXshow = min(6, n_x);
idxJ = unique(round(linspace(1, n_J, nJshow)));
idxX = unique(round(linspace(1, n_x, nXshow)));
legJ = arrayfun(@(jj) sprintf('J=%.3g', J_ref(jj)), idxJ, 'UniformOutput', false);
legX = arrayfun(@(ii) sprintf('x=%.2f', x_rel(ii)), idxX, 'UniformOutput', false);

% Re(x)
newfig('Re(x)'); hold on; grid on
for k = 1:numel(idxJ), j = idxJ(k)
    semilogy(x_rel, RE_x(:,j), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('x = r/R'); ylabel('Re'); title('Reynolds ao longo da pá')
legend(legJ,'Location','best','Box','off'); format_axes(gca)

% Re(J)
newfig('Re(J)'); hold on; grid on
for k = 1:numel(idxX), i = idxX(k)
    semilogy(J_ref, RE_x(i,:), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('J'); ylabel('Re'); title('Reynolds em estações fixas')
legend(legX,'Location','best','Box','off'); format_axes(gca)

% Mach(x)
newfig('Mach(x)'); hold on; grid on
for k = 1:numel(idxJ), j = idxJ(k)
    plot(x_rel, MA_x(:,j), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('x = r/R'); ylabel('Mach'); title('Mach ao longo da pá')
legend(legJ,'Location','best','Box','off'); format_axes(gca)

% Mach(J)
newfig('Mach(J)'); hold on; grid on
for k = 1:numel(idxX), i = idxX(k)
    plot(J_ref, MA_x(i,:), '-', 'LineWidth', 2.0, 'Color', crad(k));
end
xlabel('J'); ylabel('Mach'); title('Mach em estações fixas')
legend(legX,'Location','best','Box','off'); format_axes(gca)

% ====================== FUNÇÕES LOCAIS (utils) ==========================
function a = dcl_dalpha_local(alpha_deg, Re_u, Ma_u, Pgrid, method, iStation, jIndex, enforceMonotone)
% Inclinação local dCl/dalpha (1/deg) estimada da própria polar (dif. central).
    if nargin < 5 || isempty(method), method = 'linear'; end
    da = 0.5;  % deg
    cl_p = busca_polares_corretas(alpha_deg+da, Re_u, Ma_u, Pgrid, method, iStation, jIndex, enforceMonotone, false);
    cl_m = busca_polares_corretas(alpha_deg-da, Re_u, Ma_u, Pgrid, method, iStation, jIndex, enforceMonotone, false);
    a = (cl_p - cl_m)/(2*da);
end
