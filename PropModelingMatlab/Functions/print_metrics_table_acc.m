function print_metrics_table_acc(mCT, mCP, mETA, smape_final_vec, labels, rpm, csvname)
% Imprime e (opcional) salva CSV com N, R2, r2_shape, sMAPE acumulada.

Serie = labels(:);
N = [mCT.N; mCP.N; mETA.N];
R2 = [mCT.R2; mCP.R2; mETA.R2];
r2sh = [mCT.r2shape; mCP.r2shape; mETA.r2shape];
sMAPE = smape_final_vec(:);

T = table(Serie, N, R2, r2sh, sMAPE, ...
  'VariableNames', {'Serie', 'N', 'R2', 'r2_shape', 'sMAPE_accum_percent'});

fprintf('\n==== Comparativo Lowry × Referência – RPM %d ====\n', rpm);
disp(T);

if nargin >= 6 && ~isempty(csvname)
  try
    writetable(T, csvname);
    fprintf('-> Tabela salva em "%s"\n\n', csvname);
  catch
    warning('Não foi possível salvar CSV "%s".', csvname);
  end
end
end
