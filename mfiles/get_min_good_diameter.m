%% get max diameter
% need full psd per cruise
ff = dir('\\jett\awlab\DATA\IFCB CRUISE DATA\HOT\*\proc\ifcb_uw_psd.mat');
ff(5,:)=[]; % hot308
for i = 1:length(ff)-1 % skip 357
    load(fullfile(ff(i).folder,ff(i).name),'ifcbPsdn','ifcbPsdnCSA','diams');
   % loglog(diams,nanmax(ifcbPsdn),'k'); hold on
    [tmp indd]= max(nanmax(ifcbPsdn));
    max_good_diam(i,:) = diams(indd);

    [tmp indd]= max(nanmax(ifcbPsdnCSA));
    max_good_diam_CSA(i,:) = diams(indd);
end

save('\\jett\awlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data\max_good_diam.mat','max_good_diam','max_good_diam_CSA');