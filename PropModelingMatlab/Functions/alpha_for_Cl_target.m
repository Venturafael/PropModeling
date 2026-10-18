function alpha_deg = alpha_for_Cl_target(Cl_target, Re, Mach, polgrid)
    [NR,NM] = size(polgrid);
    ReV = zeros(NR,1); for ii=1:NR, ReV(ii) = polgrid(ii,1).Re; end
    MV  = zeros(1,NM); for jj=1:NM, MV(jj) = polgrid(1,jj).Mach; end
    [ReV,iord] = sort(ReV(:)); [MV,jord] = sort(MV(:).'); Pp = polgrid(iord,jord);

    Amins = +inf; Amaxs = -inf;
    ii = [1, max(find(ReV<=Re,1,'last')), min(find(ReV>=Re,1,'first'))]; ii = ii(~isnan(ii));
    jj = [1, max(find(MV<=Mach,1,'last')), min(find(MV>=Mach,1,'first'))]; jj = jj(~isnan(jj));
    ii = unique(ii); jj = unique(jj);
    for i = ii
        for j = jj
            if isfield(Pp(i,j),'Alpha') && ~isempty(Pp(i,j).Alpha)
                A = Pp(i,j).Alpha(:);
                Amins = min(Amins, min(A));  Amaxs = max(Amaxs, max(A));
            end
        end
    end
    if ~isfinite(Amins) || ~isfinite(Amaxs), Amins=-10; Amaxs=25; end

    f = @(a) (lookup_clarky2(a,Re,Mach,polgrid,'interp') - Cl_target).^2;
    alpha_deg = fminbnd(f, Amins, Amaxs);
end

