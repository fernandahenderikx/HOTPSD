function T = runTaxonMomentStats_fromDir(taxonFiles, vdn13, dD13, diams13, D50, D90, startBin)
% vdn13, dD13, diams13 are ALREADY restricted to startBin:end (e.g., 13:end)
% D50 and D90 are computed from vdn13 already.
% Taxon files still contain full-size vdnCSA_aloha; we truncate those to startBin:end
% before doing anything.

if nargin < 7 || isempty(startBin)
    startBin = 13;
end
epsVal = eps;

% Force shapes
diams13 = diams13(:)'; 
dD13    = dD13(:)'; 
D50     = D50(:);
D90     = D90(:);

nTax  = numel(taxonFiles);
nRows = nTax * 4;

Taxon      = strings(nRows,1);
Comparison = strings(nRows,1);
Response   = strings(nRows,1);
n          = nan(nRows,1);
Beta       = nan(nRows,1);
Beta_std   = nan(nRows,1);
R2         = nan(nRows,1);
p          = nan(nRows,1);

row = 0;

for i = 1:nTax
    
    fpath = fullfile(taxonFiles(i).folder, taxonFiles(i).name);
    S = load(fpath, 'vdnCSA_aloha');
    taxV_full = S.vdnCSA_aloha;                % [nCruise x nBins_full]
    
    % Truncate taxon PSD to match your bulk arrays (which are already 13:end)
    taxV13 = taxV_full(:, startBin:end);       % [nCruise x nBins13]

    % Truncate cruise to exclude 308, like i excluded in bulk. assume
    % cruise 308 is = 5thcruise.
    taxV13 = taxV_full([1:4 6:end], startBin:end);       % [nCruise x nBins13]
    
    % Basic sanity check: dimensions must match bulk
    if size(taxV13,2) ~= size(vdn13,2)
        error("Column mismatch after truncation for %s. taxV13 has %d cols, vdn13 has %d cols.", ...
              taxonFiles(i).name, size(taxV13,2), size(vdn13,2));
    end
    
    % Taxon label
    taxName = erase(taxonFiles(i).name, '_uw_all.mat');

    % Integrated taxon volume over truncated bins
    Vtax = trapz(diams13, taxV13, 2);          % [nCruise x 1]
    X    = log10(Vtax + epsVal);               % predictor
    
    % 1) bulk D50 vs taxon
    [b,bstd,r2,pv,nn] = fitlm_summary(X, D50);
    row=row+1; Taxon(row)=taxName; Comparison(row)="bulk"; Response(row)="D50";
    n(row)=nn; Beta(row)=b; Beta_std(row)=bstd; R2(row)=r2; p(row)=pv;
    
    % 2) bulk D90 vs taxon
    [b,bstd,r2,pv,nn] = fitlm_summary(X, D90);
    row=row+1; Taxon(row)=taxName; Comparison(row)="bulk"; Response(row)="D90";
    n(row)=nn; Beta(row)=b; Beta_std(row)=bstd; R2(row)=r2; p(row)=pv;
    
    % Remove taxon from truncated bulk and recompute moments
    bulkMinus = vdn13 - taxV13;

    Vtot = sum(bulkMinus .* dD13, 2);
    P    = (bulkMinus .* dD13) ./ (Vtot + epsVal);
    cumP = cumsum(P, 2);

    [~,i50] = min(abs(cumP - 0.5), [], 2);
    [~,i90] = min(abs(cumP - 0.9), [], 2);

    D50m = diams13(i50);
    D90m = diams13(i90);

    % 3) D50minusTax vs taxon
    [b,bstd,r2,pv,nn] = fitlm_summary(X, D50m);
    row=row+1; Taxon(row)=taxName; Comparison(row)="bulk_minus_taxon"; Response(row)="D50";
    n(row)=nn; Beta(row)=b; Beta_std(row)=bstd; R2(row)=r2; p(row)=pv;

    % 4) D90minusTax vs taxon
    [b,bstd,r2,pv,nn] = fitlm_summary(X, D90m);
    row=row+1; Taxon(row)=taxName; Comparison(row)="bulk_minus_taxon"; Response(row)="D90";
    n(row)=nn; Beta(row)=b; Beta_std(row)=bstd; R2(row)=r2; p(row)=pv;
end

T = table(Taxon, Comparison, Response, n, Beta, Beta_std, R2, p);

end


function [b, bstd, r2, pv, nn] = fitlm_summary(x, y)
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
nn = numel(y);

if nn < 3 || std(x)==0
    b = NaN; bstd = NaN; r2 = NaN; pv = NaN;
    return
end

mdl = fitlm(x, y);
b   = mdl.Coefficients.Estimate(2);
r2  = mdl.Rsquared.Ordinary;
pv  = mdl.Coefficients.pValue(2);

mdlz = fitlm(zscore(x), zscore(y));
bstd = mdlz.Coefficients.Estimate(2);
end