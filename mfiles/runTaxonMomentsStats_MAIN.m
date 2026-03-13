clear all;close all;
datadir = '\\jett\AWlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data'

% get bulk vdn by cruise
load([datadir filesep 'avg_bulk_psd_aloha_kahe.mat'],'*aloha','diams')
u=unique(cruisen_aloha);
ind = find(u==308);
u(ind)=[];
for i = 1:length(u)
    ind = find(cruisen_aloha==u(i));
    vdn(i,:) = nanmean(vdn_aloha(ind,:),1);
end
% remove anything < 7 um, remove cruise 308

vdn = vdn(:,13:end);
diams=diams(13:end);

load([datadir filesep 'binwidth.mat']);

dD = binWidth(13:end);

% normalize bulk psd
Vbulk_tot = sum(vdn .* dD, 2);
Pbulk = (vdn .* dD) ./ (Vbulk_tot + eps);

% compute D50 (median diameter)
cumP = cumsum(Pbulk, 2);
[~, idx] = min(abs(cumP - 0.5), [], 2);
D50 = diams(idx);
[~, idx] = min(abs(cumP - 0.9), [], 2);
D90 = diams(idx);


%% run single comparisons (with and without taxon removed frombulk). excluding 308.
taxonFiles = dir([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1\*_uw_all.mat']);

T = runTaxonMomentStats(taxonFiles, vdn, dD, diams, D50, D90, 13);

writetable(T, 'Taxon_Moment_Stats.csv');


%% multiple linear regressions diatoms + cyano
load([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Bacillariophyceae_uw_all.mat'], 'vdnCSA_aloha','diams')
Vdiatom_tot = trapz(diams(13:end),vdnCSA_aloha([1:4 6:end],[13:end]),2);


load([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Cyanobacteria_uw_all.mat'], 'vdnCSA_aloha','diams')
Vcyano_tot = trapz(diams(13:end),vdnCSA_aloha([1:4 6:end],[13:end]),2);



mdl_D50_vs_dia_plus_cyano = fitlm([log10(Vdiatom_tot+eps), log10(Vcyano_tot+eps)], D50);
mdl_D90_vs_dia_plus_cyano = fitlm([log10(Vdiatom_tot+eps), log10(Vcyano_tot+eps)], D90);
mdl_D50_vs_dia_plus_cyano_zscore = fitlm([zscore(log10(Vdiatom_tot+eps)), ...
             zscore(log10(Vcyano_tot+eps))], ...
             zscore(D50));
mdl_D90_vs_dia_plus_cyano_zscore = fitlm(zscore([log10(Vdiatom_tot+eps), log10(Vcyano_tot+eps)]), zscore(D90));

%% ok, multiple linear regression of all things that drive D50 and D90:
load([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Detritus_uw_all.mat'], 'vdnCSA_aloha','diams')
Vdet_tot = trapz(diams(13:end),vdnCSA_aloha([1:4 6:end],[13:end]),2);
load([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Haptophyta_uw_all.mat'], 'vdnCSA_aloha','diams')
Vhap_tot = trapz(diams(13:end),vdnCSA_aloha([1:4 6:end],[13:end]),2);

mdl_D50_vs_cyano_plus_dia_plus_det_plus_hapto = fitlm([log10(Vcyano_tot+eps), log10(Vdiatom_tot+eps), log10(Vdet_tot+eps), log10(Vhap_tot+eps)],D50);
mdl_D50_vs_cyano_plus_dia_plus_det_plus_hapto_stdz = fitlm([zscore(log10(Vcyano_tot+eps)), zscore(log10(Vdiatom_tot+eps)), zscore(log10(Vdet_tot+eps)), zscore(log10(Vhap_tot+eps))],zscore(D50));
mdl_D90_vs_cyano_plus_det = fitlm([log10(Vcyano_tot+eps), log10(Vdet_tot+eps)],D90);
mdl_D90_vs_cyano_plus_det_stdz = fitlm([zscore(log10(Vcyano_tot+eps)), zscore(log10(Vdet_tot+eps))],zscore(D90));


%% ===== Collect all mdl_* into summary tables =====
mdlNames = who('mdl_*');

ModelSummary = table();
CoefLong = table();

for k = 1:numel(mdlNames)
    name = mdlNames{k};
    mdl = eval(name); 

    % --- Model-level summary ---
    nObs = mdl.NumObservations;
    r2   = mdl.Rsquared.Ordinary;
    r2a  = mdl.Rsquared.Adjusted;

    rmse = mdl.RMSE;
    aic  = mdl.ModelCriterion.AIC;
    bic  = mdl.ModelCriterion.BIC;
    pval = mdl.ModelFitVsNullModel.Pvalue;

    ModelSummary = [ModelSummary; table( ...
        string(name), nObs, r2, r2a, rmse, aic, bic,pval, ...
        'VariableNames', {'Model','N','R2','R2_adj','RMSE','AIC','BIC','pval'})];

    % --- Coefficient-level (long-form) ---
    C = mdl.Coefficients; % table with Estimate, SE, tStat, pValue
    tmp = table();
    tmp.Model     = repmat(string(name), height(C), 1);
    tmp.Term      = string(C.Properties.RowNames);
    tmp.Estimate  = C.Estimate;
    tmp.SE        = C.SE;
    tmp.tStat     = C.tStat;
    tmp.pValue    = C.pValue;

    % handy for heatmaps:
    tmp.negLog10p = -log10(tmp.pValue + eps);

    CoefLong = [CoefLong; tmp];
end

% Sort for convenience
ModelSummary = sortrows(ModelSummary, 'R2', 'descend');
CoefLong     = sortrows(CoefLong, {'Model','pValue'}, {'ascend','ascend'});


