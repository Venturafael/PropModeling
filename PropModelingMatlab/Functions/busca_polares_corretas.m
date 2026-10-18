function [Cl,Cd,Cm,alpha_used,Re_used,Mach_used] = ...
    busca_polares_corretas(alpha_deg, Re_raw, Mach_raw, pol, method, ...
                           iStation, jIndex, enforceMonotone, doUpdate)
% Lookup robusto + "ponte" de monotonia (Re/Ma não-decrecentes por estação).
% Bilinear em (log10(Re), Mach) + PCHIP em alpha. Clamps suaves nas bordas.
%
% method: 'linear' (default) ou 'nearest'
% enforceMonotone: true => garante Re/Ma não decrescerem ao longo de J por estação
% doUpdate: true  => atualiza o histórico interno para a próxima chamada

    if nargin < 5 || isempty(method),          method='linear';    end
    if nargin < 6 || isempty(iStation),        iStation = 1;       end
    if nargin < 7 || isempty(jIndex),          jIndex = 1;         end
    if nargin < 8 || isempty(enforceMonotone), enforceMonotone = true; end
    if nargin < 9 || isempty(doUpdate),        doUpdate = true;    end
    if strcmpi(method,'interp'), method='linear'; end

    % histórico de Re/Ma usados (por estação) para garantir monotonia
    persistent RE_prev MA_prev maxStations
    if isempty(maxStations), maxStations = 0; end
    if isempty(RE_prev) || jIndex==1 || iStation>maxStations
        maxStations = max(maxStations, iStation);
        RE_prev = -inf(1, max(1024,maxStations));
        MA_prev = -inf(1, max(1024,maxStations));
    end

    Re_cand = Re_raw; Ma_cand = Mach_raw;
    if enforceMonotone
        Re_cand = max(Re_cand, RE_prev(iStation));
        Ma_cand = max(Ma_cand, MA_prev(iStation));
    end

    % --- lookup robusto ---
    [NR,NM] = size(pol);
    ReV = zeros(NR,1); for ii=1:NR, ReV(ii)=pol(ii,1).Re; end
    MV  = zeros(1,NM); for jj=1:NM, MV(jj)=pol(1,jj).Mach; end
    [ReV,iord]=sort(ReV(:)); [MV,jord]=sort(MV(:).'); P=pol(iord,jord);

    Aglob = sort(P(1,1).Alpha(:));
    a = min(max(alpha_deg, Aglob(1)+0.5), Aglob(end)-0.5);

    Re_c   = min(max(Re_cand  , ReV(1)), ReV(end));
    Mach_c = min(max(Ma_cand  , MV(1) ), MV(end));

    logReV  = log10(ReV);
    [i1,i2,t] = brk_lin(log10(Re_c), logReV);
    [j1,j2,u] = brk_lin(Mach_c     , MV);

    if strcmpi(method,'nearest')
        if t>=0.5, i1=i2; t=0; end
        if u>=0.5, j1=j2; u=0; end
    end

    [cl11,cd11,cm11] = cellAlphaPCHIP(P(i1,j1), a);
    [cl12,cd12,cm12] = cellAlphaPCHIP(P(i1,j2), a);
    [cl21,cd21,cm21] = cellAlphaPCHIP(P(i2,j1), a);
    [cl22,cd22,cm22] = cellAlphaPCHIP(P(i2,j2), a);

    % fallbacks
    if ~isfinite(cl11)||~isfinite(cd11), [cl11,cd11,cm11]=nearestNonNan(P,i1,j1,a); end
    if ~isfinite(cl12)||~isfinite(cd12), [cl12,cd12,cm12]=nearestNonNan(P,i1,j2,a); end
    if ~isfinite(cl21)||~isfinite(cd21), [cl21,cd21,cm21]=nearestNonNan(P,i2,j1,a); end
    if ~isfinite(cl22)||~isfinite(cd22), [cl22,cd22,cm22]=nearestNonNan(P,i2,j2,a); end

    % bilinear em (logRe, Mach)
    Cl = (1-t)*(1-u)*cl11 + (1-t)*u*cl12 + t*(1-u)*cl21 + t*u*cl22;
    Cd = (1-t)*(1-u)*cd11 + (1-t)*u*cd12 + t*(1-u)*cd21 + t*u*cd22;
    Cm = (1-t)*(1-u)*cm11 + (1-t)*u*cm12 + t*(1-u)*cm21 + t*u*cm22;

    alpha_used = a;
    Re_used    = 10.^((1-t)*logReV(i1) + t*logReV(i2));
    Mach_used  = (1-u)*MV(j1) + u*MV(j2);

    if doUpdate && enforceMonotone
        RE_prev(iStation) = Re_used;
        MA_prev(iStation) = Mach_used;
    end

% ---- helpers locais ----
function [i1,i2,t] = brk_lin(x,vec)
    i2 = find(vec>=x,1,'first'); if isempty(i2), i2=numel(vec); end
    if i2==1, i1=1; t=0; return, end
    i1=i2-1; v1=vec(i1); v2=vec(i2);
    t=(x-v1)/max(v2-v1,eps); t=min(max(t,0),1);
end
function [cl,cd,cm] = cellAlphaPCHIP(S, a_in)
    A=S.Alpha(:); CL=S.Cl(:); CD=S.Cd(:);
    hasCm=isfield(S,'Cm') && ~isempty(S.Cm); if hasCm, CM=S.Cm(:); end
    [A,ix]=sort(A); CL=CL(ix); CD=CD(ix); if hasCm, CM=CM(ix); end
    cl=pchip(A,CL,a_in); cd=pchip(A,CD,a_in); cm=0; if hasCm, cm=pchip(A,CM,a_in); end
end
function [cl,cd,cm]=nearestNonNan(Pmat,ii,jj,a_in)
    [I,J]=ndgrid(1:size(Pmat,1),1:size(Pmat,2)); D=hypot(I-ii,J-jj); D(ii,jj)=Inf;
    [~,k]=min(D(:)); [ii2,jj2]=ind2sub(size(D),k);
    [cl,cd,cm]=cellAlphaPCHIP(Pmat(ii2,jj2),a_in);
end
end
