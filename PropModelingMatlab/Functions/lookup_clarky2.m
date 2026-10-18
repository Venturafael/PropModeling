function [Cl,Cd,Cm,alpha_used,Re_used,Mach_used] = lookup_clarky2(alpha_deg, Re, Mach, pol, method)
% Busca CL/CD/CM a partir do banco NOVO: pol = ans.perfil (NR x NM)
% Cada célula tem: Alpha(:), Cl(:), Cd(:), Cm(:), Re (escalar), Mach (escalar)
% method: 'interp' (default) ou 'nearest'

    if nargin < 5 || isempty(method), method = 'interp'; end

    % ---- eixos Re e Mach vindos do grid ----
    [NR,NM] = size(pol);

    ReV = zeros(NR,1);
    for ii = 1:NR, ReV(ii) = pol(ii,1).Re; end

    MV  = zeros(1,NM);
    for jj = 1:NM, MV(jj) = pol(1,jj).Mach; end

    % garantir ordem crescente (caso tenha sido salvo fora de ordem)
    [ReV, iord] = sort(ReV(:));
    [MV , jord] = sort(MV(:).');
    P = pol(iord, jord);

    % limites de alpha do banco
    A0 = sort(P(1,1).Alpha(:));
    a  = min(max(alpha_deg, A0(1)), A0(end));

    % clamp Re/Mach aos limites do banco
    Re_c   = min(max(Re  , ReV(1)), ReV(end));
    Mach_c = min(max(Mach, MV(1) ), MV(end));

    % índices e frações para interp. bilinear
    [i1,i2,t] = brk_lin(Re_c  , ReV );
    [j1,j2,u] = brk_lin(Mach_c, MV  );

    if strcmpi(method,'nearest')
        if t>=0.5, i1=i2; t=0; end
        if u>=0.5, j1=j2; u=0; end
    end

    % interp. linear em alpha nas 4 células vizinhas
    [cl11,cd11,cm11] = cellAlpha(P(i1,j1), a);
    [cl12,cd12,cm12] = cellAlpha(P(i1,j2), a);
    [cl21,cd21,cm21] = cellAlpha(P(i2,j1), a);
    [cl22,cd22,cm22] = cellAlpha(P(i2,j2), a);

    % bilinear em (Re, Mach)
    Cl = (1-t)*(1-u)*cl11 + (1-t)*u*cl12 + t*(1-u)*cl21 + t*u*cl22;
    Cd = (1-t)*(1-u)*cd11 + (1-t)*u*cd12 + t*(1-u)*cd21 + t*u*cd22;
    Cm = (1-t)*(1-u)*cm11 + (1-t)*u*cm12 + t*(1-u)*cm21 + t*u*cm22;

    % retornos “usados”
    alpha_used = a;
    Re_used    = (1-t)*ReV(i1) + t*ReV(i2);
    Mach_used  = (1-u)*MV(j1)  + u*MV(j2);

    % -------- helpers --------
    function [i1,i2,t] = brk_lin(x,vec)
        i2 = find(vec>=x,1,'first'); if isempty(i2), i2=numel(vec); end
        if i2==1, i1=1; t=0; return, end
        i1 = i2-1;
        v1 = vec(i1); v2 = vec(i2);
        t  = (x - v1) / max(v2 - v1, eps);
        t  = min(max(t,0),1);
    end

    function [cl,cd,cm] = cellAlpha(S, a_in)
        A = S.Alpha(:); CL = S.Cl(:); CD = S.Cd(:);
        hasCm = isfield(S,'Cm') && ~isempty(S.Cm);
        if hasCm, CM = S.Cm(:); end

        [A,ix] = sort(A); CL=CL(ix); CD=CD(ix); if hasCm, CM=CM(ix); end
        a_loc  = min(max(a_in, A(1)), A(end));   % sem extrapolação
        cl     = interp1(A,CL,a_loc,'linear');
        cd     = interp1(A,CD,a_loc,'linear');
        cm     = 0; if hasCm, cm = interp1(A,CM,a_loc,'linear'); end
    end
end
