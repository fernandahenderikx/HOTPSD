
% cruise-avg psd and psd slopes - uncomment if need to run again. otherwise, load below.

clear all; close all;
mdir = cd;

available_files = dir(['\\jett\awlab\DATA\IFCB CRUISE DATA\HOT\HOT*\proc\ifcb_uw_psd.mat']);
psdn_aloha = [];
gmt_aloha = [];
out_all_percruise=[];

available_files(5,:)=[]; % remove hot308

for i = 1:length(available_files)%-1 % skip 356 onwards, not ready (no gps)
       load([available_files(i).folder filesep available_files(i).name],'cruise*','ifcbPsdn','ifcbTime','ifcbCount*','lat','lon','ifcbMeanDiameterBinFull*','diams');

            % grab just ALOHA:
            ind = find(lat>=22.6 & lat<=22.9 & lon>=-158.2 & lon <=-157.9);

       
        psdn_aloha(i,:) = nanmean(ifcbPsdn(ind,:),1);
        counts_aloha(i,:) = nansum(ifcbCountFull(ind,:),1);
        meandiam_aloha(i,:) = nanmean(ifcbMeanDiameterBinFull(ind,:),1);
        gmt_aloha(i,:) = nanmean(ifcbTime(ind),1);
        cruisen(i,:) = extractBetween(available_files(i).folder,'\HOT\HOT','\proc')
end




% first good IFCB diameter:
[ind, refdiam]= nanmax(psdn_aloha');
refdiam_value = diams(refdiam);
out_all=[];
for kk = 1:length(counts_aloha(:,1))
% Get log-log regression of full cruise PSD - assume weights ~ number of particles per bin
indBin   = [refdiam(kk):40]; % refdiam is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
dd = log10(psdn_aloha(kk,indBin)'); dd(dd== -Inf) = 0; % Replace -Inf with 0
mdl      = fitlm(log10(diams(indBin)')-log10(refdiam_value(kk)),dd,'Weights',counts_aloha(kk,indBin));
modelfun = @(b,x)b(1).*(x./refdiam_value(kk)).^b(2);
x1       = diams; %logspace(log10(3.1076),log10(refdiam_value(kk)),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
% plot(x1,y1,'k-')
x1       = logspace(log10(3),log10(70),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
% figure;plot(log10(meandiam_aloha(kk,indBin)'),log10(psdn_aloha(kk,indBin)'),'o','MarkerFaceColor','b')
% hold on;plot(log10(x1),log10(y1),'.')
% 
% Print out report
out = mdl.Coefficients{1:2,1:2};
% fprintf('\nSlope: %0.4f, SE: %0.4f, STD: %0.4f\nN0:10^%0.4f, SE of exponent: %0.4f, STD of exponent: %0.4f\n',...
%     -1*out(2),out(4),out(4)*sqrt(length(indBin)),out(1),out(3),out(3)*sqrt(length(indBin)))
% fprintf('\nN0(+1sigma): %i\nN0: %i\nN0(-1sigma): %i\n\n',10^(out(1)+out(3)*sqrt(length(indBin))),...
%     10^out(1),10^(out(1)-out(3)*sqrt(length(indBin))))
% 

out_all(kk,:) = [out(2,1) out(2,2)];
end

dv = datevec(gmt_aloha);
summer = find(dv(:,2)>=6&dv(:,2)<=8);
winter = find(dv(:,2)==12|dv(:,2)<=2);
spring = find(dv(:,2)>=3&dv(:,2)<=5);
fall = find(dv(:,2)>=9&dv(:,2)<=11);

disp(['summer = ' num2str(nanmean(out_all(summer,1)),2) '+/-' num2str(nanstd(out_all(summer,1)),2)])
disp(['winter = ' num2str(nanmean(out_all(winter,1)),2) '+/-' num2str(nanstd(out_all(winter,1)),2)])
disp(['spring = ' num2str(nanmean(out_all(spring,1)),2) '+/-' num2str(nanstd(out_all(spring,1)),2)])
disp(['fall = ' num2str(nanmean(out_all(fall,1)),2) '+/-' num2str(nanstd(out_all(fall,1)),2)])

nanmin(out_all(:,1))
nanmax(out_all(:,1))

sqrt(sum(out_all(summer,2).^2)/length(out_all(summer,1)))
sqrt(sum(out_all(winter,2).^2)/length(out_all(winter,1)))
sqrt(sum(out_all(fall,2).^2)/length(out_all(fall,1)))
sqrt(sum(out_all(spring,2).^2)/length(out_all(spring,1)))

