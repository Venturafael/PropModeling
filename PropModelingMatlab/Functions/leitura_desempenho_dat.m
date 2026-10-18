% =============================================================
% LEITURA_DESEMPENHO_DAT
% Curated from the TCC script of the same name. The block parser is
% unchanged. Set performance_file and target_rpm before RUN for a
% non-interactive execution. Set plot_reference_data=false to suppress
% the reader-only preview figure (model figures are not affected).
%   • Encontra o bloco de dados para um RPM qualquer (ex.: 2000,
%     8500, 20000 etc.) independentemente do número de espaços.
%   • Devolve J_exp, Pe_exp, Ct_exp, Cp_exp, target_rpm.
% =============================================================

% -------------------- arquivo de performance -----------------
if ~exist('performance_file','var') || isempty(performance_file)
  paths_local = propeller_repo_paths();
  performance_file = paths_local.performance_22x10;
end
filename = performance_file;

% ---- (opcional) mostrar RPMs existentes ---------------------
fid = fopen(filename, 'r');
assert(fid ~= -1, 'Erro ao abrir %s', filename)
allLines = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
allLines = allLines{1};

rpmTokens = regexp(allLines, 'PROP RPM\s*=\s*(\d+)', 'tokens');
rpmList = unique(cellfun(@(t) str2double(t{1}), rpmTokens(~cellfun('isempty', rpmTokens))));

fprintf('RPMs disponíveis: ');
fprintf('%d  ', rpmList);
fprintf('\n');

% -------------------- entrada do usuário ---------------------
if ~exist('target_rpm','var') || isempty(target_rpm)
  target_rpm = input('Digite o RPM desejado: ');
end

% -------------------- acha o bloco ---------------------------
expr = ['PROP RPM\s*=\s*', num2str(target_rpm)]; % regex flexível
rpm_line_idx = find(~cellfun('isempty', regexp(allLines, expr, 'once')), 1);

if isempty(rpm_line_idx)
  error('RPM %d não encontrado no arquivo.', target_rpm);
end

% -------------------- extrai linhas numéricas ----------------
data_start = rpm_line_idx + 3; % primeira linha de dados
data = [];

for k = data_start:numel(allLines)
  txt = strtrim(allLines{k});
  if isempty(txt) || ~isempty(regexp(txt, 'PROP RPM\s*=', 'once'))
    break % chegou no próximo bloco ou fim
  end
  vals = sscanf(txt, '%f');
  if numel(vals) >= 5
    data = [data; vals(2:5).']; % [J  Pe  Ct  Cp]
  end
end

if isempty(data)
  error('Bloco para RPM %d encontrado, mas sem dados.', target_rpm);
end

% -------------------- separa colunas -------------------------
J_exp = data(:, 1);
Pe_exp = data(:, 2);
Ct_exp = data(:, 3);
Cp_exp = data(:, 4);

% -------------------- plot opcional (três eixos) -------------
if ~exist('plot_reference_data','var') || isempty(plot_reference_data)
  plot_reference_data = true;
end
if plot_reference_data
figure('Color', 'w');

% eixo 1 – Ct & Cp
ax1 = axes;
plot(ax1, J_exp, Ct_exp, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Ct (Thrust)'); hold(ax1, 'on');
plot(ax1, J_exp, Cp_exp, 'k--s', 'LineWidth', 1.5, 'DisplayName', 'Cp (Power)');
xlabel(ax1, 'Avanço (J)');
ylabel(ax1, 'Ct / Cp');
grid(ax1, 'on');
ax1.YColor = 'k';
ylim(ax1, [0, max([Ct_exp; Cp_exp]) * 1.2]);
title(ax1, sprintf('Dados experimentais – RPM %d', target_rpm));

% eixo 2 – Eficiência Pe
ax2 = axes;
plot(ax2, J_exp, Pe_exp, 'r--d', 'LineWidth', 1.5, 'DisplayName', 'Pe (Eficiência)'); hold(ax2, 'on');
ax2.Color = 'none';
ax2.YAxisLocation = 'right';
ax2.XLim = ax1.XLim;
ax2.XTick = [];
ylabel(ax2, 'Eficiência (Pe)');
ylim(ax2, [0, max(Pe_exp) * 1.2]);
ax2.YColor = 'r';

% legendas
legend(ax1, {'Ct (Thrust)', 'Cp (Power)'}, 'Location', 'northwest');
legend(ax2, {'Pe (Eficiência)'}, 'Location', 'northeast');
linkaxes([ax1 ax2], 'x');
end
