% KEEP THIS
%% mean_shape_contribution_remove_group.m
% Mean-shape contribution test:
%   Compare mean normalized bulk PSD to mean normalized bulk PSD with a group removed,
%   and quantify the difference with distance metrics.
%
% Distances implemented:
%   - L1 distance: sum_d |ΔP(d)|
%   - Bray–Curtis distance: sum|A-B| / (sum(A)+sum(B))
%   - Jensen–Shannon divergence (base-2, bounded 0..1)
%
% You can run this for dinos and any other group(s) to compare mean-shape importance.
%
% REQUIRED INPUTS (in workspace):
%   Vbulk_cruise   [nCruise x nBin]   cruise-mean bulk PSD, µL L^-1 µm^-1
%   Vgroup_cruise  [nCruise x nBin]   cruise-mean PSD for the group to remove (e.g., dinos), same units
%   diams          [nBin x 1]         bin centers (µm)
%   dD             [1 x nBin]         bin widths (µm)
%
% OPTIONAL:
%   detect_limit_um (default 7)
%   remove_mode: 'subtract' (default) or 'leaveoneout_fraction' (see note)

datadir = '\\jett\AWlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data';
% get bulk vdn by cruise
load([datadir filesep 'avg_bulk_psd_aloha_kahe.mat'],'cruisen_aloha','vdn_aloha')
u=unique(cruisen_aloha);
for i = 1:length(u)
    ind = find(cruisen_aloha==u(i));
    Vbulk_cruise(i,:) = nanmean(vdn_aloha(ind,:),1);
end

%% Vgroup:
%% -------- GROUP FILES TO TEST --------
groupFiles = {
    [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Arthropoda_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Bacillariophyceae_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Chlorophyta_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Ciliophora_uw_all.mat']
      [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Cryptophyta_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Cyanobacteria_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Detritus_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Dictyochophyceae_uw_all.mat']
      [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Dinophyceae_uw_all.mat']
      [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Euglenozoa_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Haptophyta_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Radiozoa_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Unidentifiable_uw_all.mat']
     [datadir filesep 'CNN_TS4_2025-02-26_merged_Group1' filesep 'Unknown_uw_all.mat']
};

groupNames = {
    'Arthopoda'
    'Bacillariophyta'
    'Chorophyta'
    'Ciliophora'
    'Cryptophyta'
    'Cyanobacteria'
    'Detritus'
    'Dictyochophyta'
    'Dinophyta'
    'Euglenozoa'
    'Haptophyta'
    'Radiozoa'
    'Unidentifiable'    
    'Unknown'
};

nGroups = numel(groupFiles);

Vgroup_cruise = vdn_aloha;

load([datadir filesep 'binWidth']);

dD = binWidth;


%% ---------------- SETTINGS ----------------
if ~exist('detect_limit_um','var'); detect_limit_um = 7; end
remove_mode = 'subtract';   % 'subtract' is standard: bulk_minus = bulk - group

epsProb = 1e-12;            % smoothing for JSD (avoid log(0))


%% prep results
results = struct();
for g = 1:nGroups
    
    fprintf('\n==============================\n');
    fprintf('Processing group: %s\n', groupNames{g});
    
    load(groupFiles{g}, 'vdnCSA_aloha','diams','gmt_all','cruisen_all');
    
    Vgroup_cruise = vdnCSA_aloha;



%% ---------------- DEFENSIVE ALIGNMENT ----------------
diams = diams(:);
dD = dD(:)';

nBin = size(Vbulk_cruise,2);
if numel(diams) ~= nBin || numel(dD) ~= nBin
    nMin = min([nBin, numel(diams), numel(dD)]);
    warning('Aligning (trimming) to min length = %d. Check bin alignment.', nMin);
    Vbulk_cruise  = Vbulk_cruise(:,1:nMin);
    Vgroup_cruise = Vgroup_cruise(:,1:nMin);
    diams = diams(1:nMin);
    dD = dD(1:nMin);
end

% Apply cutoff
keep = diams >= detect_limit_um;
D = diams(keep);
dDk = dD(keep);

Vbulk = Vbulk_cruise(:, keep);
Vgrp  = Vgroup_cruise(:, keep);

% Sanity: ensure no negatives after subtraction (can happen from rounding/mismatch)
switch lower(remove_mode)
    case 'subtract'
        Vminus = max(Vbulk - Vgrp, 0);
    otherwise
        error('Unknown remove_mode: %s', remove_mode);
end

%% ---------------- NORMALIZE EACH CRUISE TO SHAPE (integrated fractions) ----------------
% Convert to integrated contributions and then to per-cruise fractions so rows sum to 1
Pbulk  = normalize_psd_to_fractions(Vbulk,  dDk);
Pminus = normalize_psd_to_fractions(Vminus, dDk);

% Option: restrict to cruises where both are valid distributions
ok = all(isfinite(Pbulk),2) & all(isfinite(Pminus),2) & ...
     sum(Pbulk,2)>0 & sum(Pminus,2)>0;

Pbulk_ok  = Pbulk(ok,:);
Pminus_ok = Pminus(ok,:);

fprintf('Kept %d/%d cruises after validity filtering.\n', sum(ok), size(Vbulk,1));

%% ---------------- MEAN SHAPES ----------------
Pbar_bulk  = mean(Pbulk_ok,  1, 'omitnan');
Pbar_minus = mean(Pminus_ok, 1, 'omitnan');

% Ensure they sum to 1 (numerical)
Pbar_bulk  = Pbar_bulk  / (sum(Pbar_bulk)  + eps);
Pbar_minus = Pbar_minus / (sum(Pbar_minus) + eps);


%% ---------------- DISTANCE METRICS BETWEEN MEAN SHAPES ----------------
dL1 = sum(abs(Pbar_bulk - Pbar_minus));

dBC = bray_curtis(Pbar_bulk, Pbar_minus);

dJSD = jensen_shannon_divergence(Pbar_bulk, Pbar_minus, epsProb);

fprintf('\nMean-shape distances (bulk vs bulk-without-group):\n');
fprintf('  L1 distance          = %.4f\n', dL1);
fprintf('  Bray–Curtis distance = %.4f\n', dBC);
fprintf('  Jensen–Shannon div   = %.4f (base-2)\n', dJSD);

results(g).group  = groupNames{g};
    results(g).L1     = dL1;
    results(g).BrayCurtis = dBC;
    results(g).JSD    = dJSD;
    results(g).Pbar_bulk  = Pbar_bulk;
    results(g).Pbar_minus = Pbar_minus;

end
summaryTable = struct2table(results);
disp(summaryTable(:, {'group','L1','BrayCurtis','JSD'}))

%% ---------------- OPTIONAL: PLOT MEAN SHAPES ----------------
figure('Color','w','Position',[120 120 900 360]);

subplot(1,2,1);
plot(D, Pbar_bulk,  'k-', 'LineWidth',2); hold on;
plot(D, Pbar_minus, '-', 'LineWidth',2);
set(gca,'XScale','log');
xlabel('Diameter (µm)');
ylabel('Mean normalized bulk PSD (integrated fraction per bin)');
legend({'Mean bulk','Mean bulk \ group'}, 'Location','best');
title(sprintf('Mean shapes (≥%.1f µm)', detect_limit_um));

subplot(1,2,2);
plot(D, (Pbar_minus - Pbar_bulk), 'k-', 'LineWidth',2);
set(gca,'XScale','log');
xlabel('Diameter (µm)');
ylabel('\Delta mean fraction (minus - bulk)');
yline(0,'--');
title(sprintf('\\Delta mean shape: L1=%.3f, BC=%.3f, JSD=%.3f', dL1, dBC, dJSD));

%% ========================= FUNCTIONS =========================
function P = normalize_psd_to_fractions(V, dDk)
    % V: [nCruise x nBin] in µL L^-1 µm^-1
    % dDk: [1 x nBin] in µm
    Int = V .* dDk;                 % integrated contribution per bin (µL/L)
    tot = sum(Int, 2);              % total integrated volume (µL/L)
    P = Int ./ (tot + eps);         % fractions per bin
end

function d = bray_curtis(a, b)
    a = a(:)'; b = b(:)';
    d = sum(abs(a-b)) / (sum(a) + sum(b) + eps);
end

function jsd = jensen_shannon_divergence(p, q, epsProb)
    % Base-2 Jensen–Shannon divergence, bounded 0..1 for distributions
    p = p(:)'; q = q(:)';
    p = p + epsProb; q = q + epsProb;
    p = p / sum(p);
    q = q / sum(q);
    m = 0.5*(p+q);
    jsd = 0.5*kl_div(p, m) + 0.5*kl_div(q, m);
end

function k = kl_div(p, q)
    % KL(p||q) base-2
    k = sum(p .* log2(p ./ q));
end
