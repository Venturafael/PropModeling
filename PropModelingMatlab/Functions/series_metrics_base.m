function met = series_metrics_base(y, yhat)
% Métricas entre séries:
%   - R2              : coeficiente de determinação
%   - r2shape         : Pearson^2 sobre séries centradas
%   - sMAPE_mean      : sMAPE média (%)
%   - sMAPE_acc_curve : sMAPE acumulada 0→i (%)
%   - sMAPE_acc_final : último valor da acumulada (%)

y = y(:); yhat = yhat(:);
v = isfinite(y) & isfinite(yhat);
y = y(v); yhat = yhat(v);

met.N = numel(y);
if met.N == 0
  met.R2 = NaN; met.r2shape = NaN;
  met.sMAPE_mean = NaN; met.sMAPE_acc_curve = NaN; met.sMAPE_acc_final = NaN;
  return
end

% R2
e = yhat - y;
ym = mean(y, 'omitnan');
SSres = sum(e.^2, 'omitnan');
SStot = sum((y - ym).^2, 'omitnan');
met.R2 = 1 - SSres/max(SStot, eps);

% r2_shape (Pearson^2 centrado)
met.r2shape = local_pearson_r2(yhat, y);

% sMAPE média e acumulada
term = 2*abs(e) ./ (abs(yhat) + abs(y) + eps); % termos em fração
met.sMAPE_mean = 100*mean(term, 'omitnan');
met.sMAPE_acc_curve = 100*cumsum(term)./(1:met.N).';
met.sMAPE_acc_final = met.sMAPE_acc_curve(end);
end

% ---- helper local
function r2 = local_pearson_r2(x, y)
v = isfinite(x) & isfinite(y); x = x(v); y = y(v);
if numel(x) < 2, r2 = NaN; return, end
x = x - mean(x); y = y - mean(y);
sx2 = sum(x.^2); sy2 = sum(y.^2);
if sx2 <= 0 || sy2 <= 0, r2 = NaN; return, end
r = (x'*y)/sqrt(sx2*sy2);
r2 = r.^2;
end
