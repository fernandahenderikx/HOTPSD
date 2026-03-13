% KEEP FOR PAPER
%Is there seasonal variability in the PSD? (overall test)
% dinos:
datadir = '\\jett\AWlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data'

files = dir([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1\*_uw_all.mat'])
TABLE={};
for mmm = 1:length(files)
    load(fullfile(files(mmm).folder, files(mmm).name),'vdnCSA_aloha','diams','gmt_all','cruisen_all')
    ind = find(cruisen_all==308);
    vdnCSA_aloha(ind,:)=[];
    gmt_all(ind)=[];
    % make season array
    if ~isa(gmt_all,'datetime')
        gmt_all = datetime(gmt_all,'ConvertFrom','datenum');
    end
    mo = month(gmt_all);
    season = strings(size(mo));
    season(ismember(mo,[12 1 2])) = "Winter";
    season(ismember(mo,[3 4 5]))  = "Spring";
    season(ismember(mo,[6 7 8]))  = "Summer";
    season(ismember(mo,[9 10 11]))= "Fall";
    season = categorical(season, ["Winter","Spring","Summer","Fall"]);
    
    V = vdnCSA_aloha;
    % remove all data < 7 um
    keep7 = (diams >= 7);
    V7 = V(:, keep7);
    diams7 = diams(keep7);
     % Bin widths (µm) - from your PSD config
    load([datadir filesep 'binwidth.mat']);
    dD = binWidth(:)'; % 1 x nBin_expected
    dD7 = dD(keep7);

    % -------- Magnitude metric: integrated volume >= detect_limit_um
    Vtot7 = sum(V7 .* dD7, 2);                 % integrated volume (µL/L)
    bulk = log10(Vtot7 + eps);                 % log-transform
    
    %One-way ANOVA + Tukey post-hoc (MATLAB)
    [p, tbl, stats] = anova1(bulk, season);
    [c, m] = multcompare(stats);

    sig = c(c(:,end) < 0.05, :);
    parts = strings(size(sig,1),1);
    for k=1:size(sig,1)
        i = sig(k,1); j = sig(k,2);
        diff = sig(k,4);
        p = sig(k,6);
        if diff > 0
            parts(k) = sprintf('%s > %s (p=%.3g)', stats.gnames{i}, stats.gnames{j}, p);
        else
            parts(k) = sprintf('%s > %s (p=%.3g)', stats.gnames{j}, stats.gnames{i}, p);
        end
    end
    txt = strjoin(parts, '; ');

    res= {files(mmm).name tbl{2,5} tbl{2,3} tbl{3,3} p txt};

    TABLE=[TABLE; res];


end



