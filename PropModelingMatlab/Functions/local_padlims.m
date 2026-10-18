function L = local_padlims(v, lo, hi)
% Expande limites [min max] com frações lo/hi do range.

v = v(isfinite(v));
if isempty(v), L = [0 1]; return, end
vr = max(v) - min(v);
if vr <= 0
  L = [min(v) - eps, max(v) + eps];
else
  L = [min(v) - lo*vr, max(v) + hi*vr];
end
end
