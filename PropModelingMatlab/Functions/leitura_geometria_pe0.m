% LEITURA_GEOMETRIA_PE0  Lê geometria APC .PE0 para o workspace chamador.
%
% Curated from the TCC script of the same name. The numerical parser is
% unchanged; only the hard-coded path and destructive CLEARVARS were
% removed. Set geometry_file before RUN to choose another propeller.
% Fluxo: abrir arquivo -> coletar linhas -> localizar bloco numérico ->
% ler 13 colunas -> separar vetores -> extrair RADIUS/BLADES.

if ~exist('geometry_file','var') || isempty(geometry_file)
  paths_local = propeller_repo_paths();
  geometry_file = paths_local.geometry_22x10;
end
filename = geometry_file;

fid = fopen(filename, 'r');
if fid == -1
  error('Erro ao abrir o arquivo: %s', filename);
end

lines = {};
while ~feof(fid)
  lines{end+1} = fgetl(fid); %#ok<AGROW>
end
fclose(fid);

% Localiza o bloco de dados (entre cabeçalho das colunas e a linha com RADIUS).
hdr_idx = find(contains(lines, 'STATION') & contains(lines, 'CHORD'), 1, 'first');
radius_idx = find(contains(lines, 'RADIUS:'), 1, 'first');

if isempty(hdr_idx) || isempty(radius_idx)
  error('Cabeçalho de colunas ou linha de RADIUS não encontrados.');
end

start_idx = hdr_idx + 2;
end_idx = radius_idx - 1;

% Lê o bloco numérico. Aceita somente linhas com 13 valores.
raw_data = lines(start_idx:end_idx);
data = zeros(numel(raw_data), 13);
k = 0;

for i = 1:numel(raw_data)
  s = strtrim(raw_data{i});
  if isempty(s)
    continue
  end
  % Normaliza vírgula decimal, se houver.
  s = strrep(s, ',', '.');

  row = sscanf(s, '%f');
  if numel(row) == 13
    k = k + 1;
    data(k, :) = row.';
  end
end

if k == 0
  error('Nenhuma linha numérica válida (13 colunas) foi encontrada em %s.', filename);
end

data = data(1:k, :);

% Separa as colunas.
STATION = data(:, 1);
CHORD = data(:, 2);
PITCH_QUOTED = data(:, 3);
PITCH_LE_TE = data(:, 4);
PITCH_PRATHER = data(:, 5);
SWEEP = data(:, 6);
THICKNESS_RATIO = data(:, 7);
TWIST = data(:, 8);
MAX_THICK = data(:, 9);
CROSS_SECTION = data(:, 10);
ZHIGH = data(:, 11);
CGY = data(:, 12);
CGZ = data(:, 13);

% Extrai RADIUS e BLADES do cabeçalho.
line_radius = lines{radius_idx};
line_blades = lines{find(contains(lines, 'BLADES:'), 1, 'first')};

% Extração com tolerância a notação científica.
tokR = regexp(line_radius, 'RADIUS:\s*([0-9Ee+\-\.]+)', 'tokens', 'once');
tokB = regexp(line_blades, 'BLADES:\s*(\d+)', 'tokens', 'once');

RADIUS = NaN;
BLADES = NaN;
if ~isempty(tokR), RADIUS = str2double(tokR{1}); end
if ~isempty(tokB), BLADES = str2double(tokB{1}); end
