function v = linfill1d(v)
% Preenche NaNs por interpolação 1D (com extrapolação às bordas).
    if isrow(v), v=v(:); end
    idx = find(isfinite(v));
    if isempty(idx), v=zeros(size(v)); return, end
    x=(1:numel(v))'; vv=v;
    vv(isnan(vv)) = interp1(x(idx), v(idx), x(isnan(vv)), 'linear', 'extrap');
    v=vv;
end
