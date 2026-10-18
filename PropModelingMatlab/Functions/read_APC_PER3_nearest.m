function E = read_APC_PER3_nearest(filename, rpm_pref)
    fid = fopen(filename,'r');  assert(fid~=-1,'Erro ao abrir %s',filename);
    LL = textscan(fid,'%s','Delimiter','\n'); fclose(fid); L = LL{1};
    rpmTokens = regexp(L,'PROP RPM\s*=\s*(\d+)','tokens');
    rpmList   = unique( cellfun(@(t) str2double(t{1}), rpmTokens(~cellfun('isempty',rpmTokens))) );
    fprintf('RPMs disponíveis:  '); fprintf('%d  ',rpmList); fprintf('\n');
    [~,ix] = min(abs(rpmList - rpm_pref));
    target_rpm = rpmList(ix);
    fprintf('Usando RPM experimental mais próximo: %d (pref=%d)\n', target_rpm, rpm_pref);
    expr = ['PROP RPM\s*=\s*' , num2str(target_rpm)];
    rpm_line_idx = find(~cellfun('isempty', regexp(L, expr)), 1);
    assert(~isempty(rpm_line_idx),'Bloco do RPM não encontrado.');
    data_start = rpm_line_idx + 3;
    data = [];
    for k = data_start:numel(L)
        txt = strtrim(L{k});
        if isempty(txt) || ~isempty(regexp(txt,'PROP RPM\s*=', 'once')), break; end
        vals = sscanf(txt,'%f');
        if numel(vals) >= 5
            data = [data ; vals(2:5).'];  % [J  Pe  Ct  Cp]
        end
    end
    assert(~isempty(data),'Bloco encontrado, mas sem dados.');
    E.J  = data(:,1);
    E.Pe = data(:,2);
    E.Ct = data(:,3);
    E.Cp = data(:,4);
    E.target_rpm = target_rpm;
end
