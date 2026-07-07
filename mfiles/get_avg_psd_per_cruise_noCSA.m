%% plot HOT PSDs
clear all; close all;
mdir = cd;

available_files = dir(['\\jett\awlab\DATA\IFCB CRUISE DATA\HOT\HOT*\proc\ifcb_uw_psd.mat']);

psdn_aloha = [];
vdn_aloha = [];
psdn_kahe = [];
vdn_kahe = [];
gmt_aloha = [];
gmt_kahe = [];
cruisen_aloha = [];
cruisen_kahe = [];
meddiam_aloha=[];
meddiam_kahe=[];
for i = 1:length(available_files) %
       load([available_files(i).folder filesep available_files(i).name],'ifcbPsdn*','ifcbTime','ifcbVDn*','lat','lon','diams','ifcbMedDiam');
    c = str2num(char(extractBetween([available_files(i).folder filesep available_files(i).name],'\HOT\HOT','\proc')));
            % grab just ALOHA:
            ind = find(lat>=22.6 & lat<=22.9 & lon>=-158.2 & lon <=-157.9);

            % grab just coastal/KAHE:
            ind2 = find(lat>=21.3& lat<=21.8 & lon>=-158.5 & lon <=-158.2);

        
        psdn_aloha = [psdn_aloha; ifcbPsdn(ind,:)];
        vdn_aloha = [vdn_aloha; ifcbVDn(ind,:)];

        psdn_kahe = [psdn_kahe; ifcbPsdn(ind2,:)];
        vdn_kahe = [vdn_kahe; ifcbVDn(ind2,:)];

        gmt_aloha = [gmt_aloha; ifcbTime(ind)];
        gmt_kahe = [gmt_kahe; ifcbTime(ind2)];

        cruisen_aloha = [cruisen_aloha; repmat(c,length(ind),1)];
        cruisen_kahe = [cruisen_kahe; repmat(c,length(ind2),1)];

end
psdn_aloha_composite = nanmean(psdn_aloha);
psdn_kahe_composite = nanmean(psdn_kahe);
vdn_aloha_composite = nanmean(vdn_aloha);
vdn_kahe_composite = nanmean(vdn_kahe);

% aloha summer and winter
dv = datevec(gmt_aloha);
summer = find(dv(:,2)>=6&dv(:,2)<=8);
winter = find(dv(:,2)>=12|dv(:,2)<=2);
fall = find(dv(:,2)>=9&dv(:,2)<=11);
spring = find(dv(:,2)>=3&dv(:,2)<=5);
allbutsummer = find(dv(:,2)<=5|dv(:,2)>=9);

psdn_aloha_summer = psdn_aloha(summer,:);
psdn_aloha_winter = psdn_aloha(winter,:);
psdn_aloha_fall = psdn_aloha(fall,:);
psdn_aloha_spring = psdn_aloha(spring,:);
psdn_aloha_allbutsummer = psdn_aloha(allbutsummer,:);
vdn_aloha_summer = vdn_aloha(summer,:);
vdn_aloha_winter = vdn_aloha(winter,:);
vdn_aloha_fall = vdn_aloha(fall,:);
vdn_aloha_spring = vdn_aloha(spring,:);
vdn_aloha_allbutsummer = vdn_aloha(allbutsummer,:);
gmt_aloha_summer = gmt_aloha(summer);
gmt_aloha_winter = gmt_aloha(winter);
gmt_aloha_fall = gmt_aloha(fall);
gmt_aloha_spring = gmt_aloha(spring);
gmt_aloha_allbutsummer = gmt_aloha(allbutsummer);
cruisen_aloha_summer = cruisen_aloha(summer);
cruisen_aloha_winter = cruisen_aloha(winter);
cruisen_aloha_fall = cruisen_aloha(fall);
cruisen_aloha_spring = cruisen_aloha(spring);
cruisen_aloha_allbutsummer = cruisen_aloha(allbutsummer);


% kahe summer and winter
dv = datevec(gmt_kahe);
summer = find(dv(:,2)>=6&dv(:,2)<=8);
winter = find(dv(:,2)>=12|dv(:,2)<=2);
psdn_kahe_summer = psdn_kahe(summer,:);
psdn_kahe_winter = psdn_kahe(winter,:);
vdn_kahe_summer = vdn_kahe(summer,:);
vdn_kahe_winter = vdn_kahe(winter,:);
psdn_kahe_summer_composite = nanmean(psdn_kahe(summer,:));
psdn_kahe_winter_composite = nanmean(psdn_kahe(winter,:));
vdn_kahe_summer_composite = nanmean(vdn_kahe(summer,:));
vdn_kahe_winter_composite = nanmean(vdn_kahe(winter,:));
gmt_kahe_summer = gmt_kahe(summer);
gmt_kahe_winter = gmt_kahe(winter);
cruisen_kahe_summer = cruisen_kahe(summer);
cruisen_kahe_winter = cruisen_kahe(winter);


info  = {'summer = find(dv(:,2)>=6&dv(:,2)<=8); winter = find(dv(:,2)>=12|dv(:,2)<=2)'};
save('\\jett\awlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data\avg_bulk_psd_aloha_kahe_notCSA.mat','psdn_aloha*','vdn_aloha*','psdn_kahe*','vdn_kahe*', 'info','diams','gmt*','cruisen*')

