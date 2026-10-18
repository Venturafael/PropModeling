function S = read_APC_PE0(filename)
    fid = fopen(filename,'r');
    if fid==-1, error('Erro ao abrir %s',filename); end
    C = textscan(fid,'%s','Delimiter','\n'); fclose(fid); L = C{1};
    start_idx = find(contains(L,'STATION')&contains(L,'CHORD')) + 2;
    end_idx   = find(contains(L,'RADIUS:'),1) - 1;
    assert(~isempty(start_idx)&&~isempty(end_idx),'Bloco numérico não encontrado em %s',filename);
    raw = L(start_idx:end_idx);
    data = zeros(numel(raw),13); k=0;
    for i=1:numel(raw)
        row = sscanf(raw{i},'%f');
        if numel(row)==13, k=k+1; data(k,:)=row.'; end
    end
    data = data(1:k,:);
    S.STATION         = data(:,1);
    S.CHORD           = data(:,2);
    S.PITCH_QUOTED    = data(:,3);
    S.PITCH_LE_TE     = data(:,4);
    S.PITCH_PRATHER   = data(:,5);
    S.SWEEP           = data(:,6);
    S.THICKNESS_RATIO = data(:,7);
    S.TWIST           = data(:,8);
    S.MAX_THICK       = data(:,9);
    S.CROSS_SECTION   = data(:,10);
    S.ZHIGH           = data(:,11);
    S.CGY             = data(:,12);
    S.CGZ             = data(:,13);
    line_radius = L{find(contains(L,'RADIUS:'),1)};
    line_blades = L{find(contains(L,'BLADES:'),1)};
    S.RADIUS = str2double(regexp(line_radius,'RADIUS:\s*([\d.]+)','tokens','once'));
    S.BLADES = str2double(regexp(line_blades,'BLADES:\s*(\d+)','tokens','once'));
end

