%% Paper Figures and Tables

clear all; close all;
figsdir = '\\jett\awlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\figures';
datadir = '\\jett\AWlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data'
load([datadir filesep 'nsamples_percruise.mat'])
load([datadir filesep 'max_good_diam.mat'])

%% prep data

clearvars -except max_good_diam*  good_cruises figsdir datadir
f = dir([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1'])
f([1:2 13 15])=[]; % exclude "bad" data and forams (no data)
%f([1:2])=[];
% edit some labels
Labels = extractBefore({f.name},'_uw_all.mat'); Labels{1,8}='Dictyochophyta'; Labels{1,9}='Dinoflagellata'; Labels{1,2}='Bacillariophyta';


good_cruises = [1:4 6:45]; % remove 308, only 2 samples

for i = 1:length(f)
    load(fullfile(f(i).folder,f(i).name),'vdn_aloha','psdn_aloha','diams','gmt_all');
    
    vdnCSA_aloha=vdn_aloha(good_cruises,:);
    psdCSA_aloha=psdn_aloha(good_cruises,:);
    gmt_all=gmt_all(good_cruises);


    % force something here: depending on whether cruise time comes from the
    % minimum time or average time, cruise HOT339 becomes either summer or
    % fall. force it to be August since sampling was at the end of August.
    dv=datevec(gmt_all-10/24);
    %summer = find(dv(:,2)>=6&dv(:,2)<=8);
    winter = find(dv(:,2)==12|dv(:,2)<=2);
    spring = find(dv(:,2)>=3&dv(:,2)<=5);
    %fall = find(dv(:,2)>=9&dv(:,2)<=11);

    summer =[1
     8
     9
    10
    15
    16
    22
    25 % added this one, because it started in August but finished in sep
    27
    34
    35
];

    fall = [3
     4
    11
    17
    18
    28
    29
    30
    31
    36
    37
    41
    42
    43
    44];


    ALL_csa_summer(i,:) = nanmean(vdnCSA_aloha(summer,:));
    ALL_csa_winter(i,:) = nanmean(vdnCSA_aloha(winter,:));
    ALL_csa_spring(i,:) = nanmean(vdnCSA_aloha(spring,:));
    ALL_csa_fall(i,:) = nanmean(vdnCSA_aloha(fall,:));
    ALL_csa_summer_std(i,:) = nanstd(vdnCSA_aloha(summer,:))./sqrt(length(summer));
    ALL_csa_winter_std(i,:) = nanstd(vdnCSA_aloha(winter,:))./sqrt(length(winter));
    ALL_csa_spring_std(i,:) = nanstd(vdnCSA_aloha(spring,:))./sqrt(length(spring));
    ALL_csa_fall_std(i,:) = nanstd(vdnCSA_aloha(fall,:))./sqrt(length(fall));

    ALL_volconc(i,:) =  trapz(diams,vdnCSA_aloha,2);

    abund_summer(i,:) = nanmean(trapz(diams,psdCSA_aloha(summer,:),2));
    abund_summer_se(i,:) = nanstd(trapz(diams,psdCSA_aloha(summer,:),2))./sqrt(length(summer));
    abund_winter(i,:) = nanmean(trapz(diams,psdCSA_aloha(winter,:),2));
    abund_winter_se(i,:) = nanstd(trapz(diams,psdCSA_aloha(winter,:),2))./sqrt(length(winter));
    abund_spring(i,:) = nanmean(trapz(diams,psdCSA_aloha(spring,:),2));
    abund_spring_se(i,:) = nanstd(trapz(diams,psdCSA_aloha(spring,:),2))./sqrt(length(spring));    
    abund_fall(i,:) = nanmean(trapz(diams,psdCSA_aloha(fall,:),2));
    abund_fall_se(i,:) = nanstd(trapz(diams,psdCSA_aloha(fall,:),2))./sqrt(length(fall));
end

f2 = dir([datadir filesep 'CNN_TS4_2025-02-26_merged_Group3'])
f2([1:2 4 37])=[]; % no forams, meta, etc


for i = 1:length(f2)
    load(fullfile(f2(i).folder,f2(i).name),'vdn_aloha','diams','gmt_all');
    %    vdnCSA_aloha(vdnCSA_aloha==0)=NaN;

    vdnCSA_aloha=vdn_aloha(good_cruises,:);
    gmt_all=gmt_all(good_cruises);

    ALL_csa_g2(i,:) = nanmean(vdnCSA_aloha)  ;

    dv=datevec(gmt_all-10/24);
    %summer = find(dv(:,2)>=6&dv(:,2)<=8);
    winter = find(dv(:,2)==12|dv(:,2)<=2);
    spring = find(dv(:,2)>=3&dv(:,2)<=5);
    %fall = find(dv(:,2)>=9&dv(:,2)<=11);
 summer =[1
     8
     9
    10
    15
    16
    22
    25 % added this one
    27
    34
    35
];

    fall = [3
     4
    11
    17
    18
    28
    29
    30
    31
    36
    37
    41
    42
    43
    44];


    ALL_csa_summer_g2(i,:) = nanmean(vdnCSA_aloha(summer,:));
    ALL_csa_winter_g2(i,:) = nanmean(vdnCSA_aloha(winter,:));
     ALL_csa_spring_g2(i,:) = nanmean(vdnCSA_aloha(spring,:));
    ALL_csa_fall_g2(i,:) = nanmean(vdnCSA_aloha(fall,:));
    ALL_csa_summer_std_g2(i,:) = nanstd(vdnCSA_aloha(summer,:))./sqrt(length(summer));
    ALL_csa_winter_std_g2(i,:) = nanstd(vdnCSA_aloha(winter,:))./sqrt(length(winter));
    ALL_csa_spring_std_g2(i,:) = nanstd(vdnCSA_aloha(spring,:))./sqrt(length(spring));
    ALL_csa_fall_std_g2(i,:) = nanstd(vdnCSA_aloha(fall,:))./sqrt(length(fall));

end

tot_volconc = nansum(ALL_volconc);




%% Fig. 1: Cruise-averaged particle size distribution 
load([datadir filesep 'avg_bulk_psd_aloha_kahe_notCSA.mat'],'psdn_aloha','vdn_aloha','gmt_aloha','cruisen_aloha')

u=unique(cruisen_aloha);
ind = find(u==308);
u(ind)=[];
clear XX t
for i = 1:length(u)
    ind = find(cruisen_aloha==u(i))
    XX(i,:)=nanmean(vdn_aloha(ind,:),1); %X(X==0)=NaN;
    yy(i,:) = nanmean(psdn_aloha(ind,:),1);
    t(i,:) = datevec(nanmin(gmt_aloha(ind))); % min or mean will give different results (some cruises started in summer, ended in fall)
end

% NOW do seasons. these will now match the hard coded ones above.
summer = find(t(:,2)>=6&t(:,2)<=8); %11
winter = find(t(:,2)==12|t(:,2)<=2); % 11
spring = find(t(:,2)>=3&t(:,2)<=5); %7
fall = find(t(:,2)>=9&t(:,2)<=11);%15

figure
subplot 122

X=XX(spring,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xm=Xm.*diams.*log(10); % to allow area comparison.
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    Xs=Xs.*diams.*log(10); % to allow area comparison.
    h1=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
X=XX(summer,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xm=Xm.*diams.*log(10); % to allow area comparison.
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    Xs=Xs.*diams.*log(10); % to allow area comparison.
    h2=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('green'),'LineWidth',2},0.5); hold on
X=XX(fall,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xm=Xm.*diams.*log(10); % to allow area comparison.
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    Xs=Xs.*diams.*log(10); % to allow area comparison.
    h3=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('orange'),'LineWidth',2},0.5); hold on
X=XX(winter,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xm=Xm.*diams.*log(10); % to allow area comparison.
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    Xs=Xs.*diams.*log(10); % to allow area comparison.
    h4=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('black'),'LineWidth',2},0.5); hold on
    ylim([0 0.11])
xlim([4 100])
set(gca,'XScale','log','YScale','linear')

 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
xlabel('Diameter, \mum')
ylabel('\muL L^-^1 per bin')
set(gca,'FontSize',16)
addlabel('b)')


subplot 121
% first plot all grey lines for each cruise
plot(diams,yy,'Color',[.7 .7 .7]); hold on

X=yy(spring,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    h1=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
X=yy(summer,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    h2=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('green'),'LineWidth',2},0.5); hold on
X=yy(fall,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    h3=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('orange'),'LineWidth',2},0.5); hold on
X=yy(winter,:);
    Xm=nanmean(X,1);% X(X==0)=NaN;
    Xs=nanstd(X,1); %Xs(Xs==0)=NaN;
    h4=shadedErrorBar(diams,Xm,Xs,{'Color',rgb('black'),'LineWidth',2},0.5); hold on

set(gca,'XScale','log','YScale','log')
ylim([10^-1 10^5])
xlim([4 100])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
xlabel('Diameter, \mum')
ylabel('# L^-^1 \mum^-^1')
set(gca,'FontSize',16)
addlabel('a)')
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],'Spring','Summer','Fall','Winter')

 set(gcf,'Position',[413         391        1253         531]);

savefig(gcf,[figsdir filesep 'Fig1.fig']); 
print(gcf,[figsdir filesep 'Fig1'],'-dpng','-r300'); 



%% Fig 2: average volume distributions partitioned by major group


figure(22222)
colors = cmocean('balance',length(f));
subplot 141
X = ALL_csa_spring.*diams.*log(10);
h=area(diams,X');
for j = 1:length(h)
h(j).FaceColor = colors(j,:);   %
end
%ylim([0 100])
hold on
set(gca,'XScale','log')
ylim([0 0.06])
xlim([5 100])
 vline(nanmean(max_good_diam),'--K')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 title('Spring')
 ylabel('\muL L^-^1 per bin')
 xlabel('Diameter, \mum')
addlabel('a)')

subplot 142
X = ALL_csa_summer.*diams.*log(10);
h=area(diams,X');
for j = 1:length(h)
h(j).FaceColor = colors(j,:);   %
end
%ylim([0 100])
hold on
set(gca,'XScale','log')
ylim([0 0.06])
xlim([5 100])
 vline(nanmean(max_good_diam),'--K')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 title('Summer')
% ylabel('\muL L^-^1 per bin')
set(gca,'YTickLabel',[])
 xlabel('Diameter, \mum')
 addlabel('b)')

subplot 143
X = ALL_csa_fall.*diams.*log(10);
h=area(diams,X');
for j = 1:length(h)
h(j).FaceColor = colors(j,:);   %
end
%ylim([0 100])
hold on
set(gca,'XScale','log')
ylim([0 0.06])
xlim([5 100])
 vline(nanmean(max_good_diam),'--K')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 title('Fall')
 %ylabel('\muL L^-^1 per bin')
 set(gca,'YTickLabel',[])
 xlabel('Diameter, \mum')
 addlabel('c)')

subplot 144
X = ALL_csa_winter.*diams.*log(10);
h=area(diams,X');
for j = 1:length(h)
h(j).FaceColor = colors(j,:);   %
end
%ylim([0 100])
hold on
set(gca,'XScale','log')
ylim([0 0.06])
xlim([5 100])
 vline(nanmean(max_good_diam),'--K')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 title('Winter')
 %ylabel('\muL L^-^1 per bin')
 set(gca,'YTickLabel',[])
 xlabel('Diameter, \mum')
addlabel('d)')
set(gcf,'Position',[198         444        1599         388])
%leg=legend(Labels,'Orientation', 'horizontal', 'Location', 'SouthOutside'); leg.NumColumns = 3;
%savefig(gcf,[figsdir filesep 'Fig2_wLeg.fig']); 
%print(gcf,[figsdir filesep 'Fig2_wLeg'],'-dpng','-r300'); 
savefig(gcf,[figsdir filesep 'Fig2.fig']); 
print(gcf,[figsdir filesep 'Fig2'],'-dpng','-r300'); 




 %% Fig 3 - average volume concentration per season portioned by major taxonomic group 


figure(22211)
data = [trapz(diams,ALL_csa_spring,2)'; ...
        trapz(diams,ALL_csa_summer,2)'; ...
        trapz(diams,ALL_csa_fall,2)'; ...
        trapz(diams,ALL_csa_winter,2)'];
data_plot = fliplr(data);   % <-- invert stacking order
h = bar(data_plot,'stacked','FaceColor','flat');
colors = cmocean('balance',length(ALL_csa_spring(:,1)));

for j = 1:length(h)
    h(j).CData = colors(end-j+1,:);   % <-- reverse colors to keep taxa colors unchanged
end

set(gca,'XTick',[1 2 3 4], ...
        'XTickLabel',{'spring','summer','fall','winter'}, ...
        'FontSize',16)

xlim([0.5 4.5])
ylim([0 0.045])
ylabel('volume concentration, \muL L^-1')
text(0.02,0.98,'a)','Units','Normalized', ...
    'VerticalAlignment','Top','FontSize',16)
leg = legend(h(end:-1:1), Labels, 'Orientation','horizontal','Location','east'); leg.NumColumns = 1;
savefig(gcf,[figsdir filesep 'Fig3_wLeg.fig']); 
print(gcf,[figsdir filesep 'Fig3_wLeg'],'-dpng','-r300'); % after manually fixing things




%% Fig. 4: Average volumetric size distributions for 6 major taxonomic groups per season

figure
subplot(2,3,1)
U=2;
shadedErrorBar(diams,ALL_csa_spring(U,:).*diams.*log(10),ALL_csa_spring_std(U,:).*diams.*log(10),{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
shadedErrorBar(diams,ALL_csa_summer(U,:).*diams.*log(10),ALL_csa_summer_std(U,:).*diams.*log(10),{'Color',rgb('green'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_fall(U,:).*diams.*log(10),ALL_csa_fall_std(U,:).*diams.*log(10),{'Color',rgb('orange'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_winter(U,:).*diams.*log(10),ALL_csa_winter_std(U,:).*diams.*log(10),{'Color',rgb('black'),'LineWidth',2},0.5); 
set(gca,'XScale','log'); axis tight; xlim([3 100])
xlabel('diamemeter, \mum'); ylabel('\muL L^-^1 per bin')
title('Bacillariophyta')
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 addlabel('a)')
 subplot(2,3,2)
U=9;
shadedErrorBar(diams,ALL_csa_spring(U,:).*diams.*log(10),ALL_csa_spring_std(U,:).*diams.*log(10),{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
shadedErrorBar(diams,ALL_csa_summer(U,:).*diams.*log(10),ALL_csa_summer_std(U,:).*diams.*log(10),{'Color',rgb('green'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_fall(U,:).*diams.*log(10),ALL_csa_fall_std(U,:).*diams.*log(10),{'Color',rgb('orange'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_winter(U,:).*diams.*log(10),ALL_csa_winter_std(U,:).*diams.*log(10),{'Color',rgb('black'),'LineWidth',2},0.5); 
set(gca,'XScale','log'); axis tight; xlim([3 100])
xlabel('diamemeter, \mum');ylabel('\muL L^-^1 per bin')
title('Dinoflagellata')
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 addlabel('b)')
 subplot(2,3,3)
U=6;
shadedErrorBar(diams,ALL_csa_spring(U,:).*diams.*log(10),ALL_csa_spring_std(U,:).*diams.*log(10),{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
shadedErrorBar(diams,ALL_csa_summer(U,:).*diams.*log(10),ALL_csa_summer_std(U,:).*diams.*log(10),{'Color',rgb('green'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_fall(U,:).*diams.*log(10),ALL_csa_fall_std(U,:).*diams.*log(10),{'Color',rgb('orange'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_winter(U,:).*diams.*log(10),ALL_csa_winter_std(U,:).*diams.*log(10),{'Color',rgb('black'),'LineWidth',2},0.5); 
set(gca,'XScale','log'); axis tight; xlim([3 100])
xlabel('diamemeter, \mum'); ylabel('\muL L^-^1 per bin')
title('Cyanobacteria')
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 addlabel('c)')
 subplot(2,3,4)
U=7;
shadedErrorBar(diams,ALL_csa_spring(U,:).*diams.*log(10),ALL_csa_spring_std(U,:).*diams.*log(10),{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
shadedErrorBar(diams,ALL_csa_summer(U,:).*diams.*log(10),ALL_csa_summer_std(U,:).*diams.*log(10),{'Color',rgb('green'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_fall(U,:).*diams.*log(10),ALL_csa_fall_std(U,:).*diams.*log(10),{'Color',rgb('orange'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_winter(U,:).*diams.*log(10),ALL_csa_winter_std(U,:).*diams.*log(10),{'Color',rgb('black'),'LineWidth',2},0.5); 
set(gca,'XScale','log'); axis tight; xlim([3 100])
xlabel('diamemeter, \mum'); ylabel('\muL L^-^1 per bin')
title('Detritus')
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 addlabel('d)')
 subplot(2,3,5)
U=13;
shadedErrorBar(diams,ALL_csa_spring(U,:).*diams.*log(10),ALL_csa_spring_std(U,:).*diams.*log(10),{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
shadedErrorBar(diams,ALL_csa_summer(U,:).*diams.*log(10),ALL_csa_summer_std(U,:).*diams.*log(10),{'Color',rgb('green'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_fall(U,:).*diams.*log(10),ALL_csa_fall_std(U,:).*diams.*log(10),{'Color',rgb('orange'),'LineWidth',2},0.5); 
shadedErrorBar(diams,ALL_csa_winter(U,:).*diams.*log(10),ALL_csa_winter_std(U,:).*diams.*log(10),{'Color',rgb('black'),'LineWidth',2},0.5); 
set(gca,'XScale','log'); axis tight; xlim([3 100])
xlabel('diamemeter, \mum'); ylabel('\muL L^-^1 per bin')
title('Unidentifiable')
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 addlabel('e)')
 subplot(2,3,6)
U=11;
h1=shadedErrorBar(diams,ALL_csa_spring(U,:).*diams.*log(10),ALL_csa_spring_std(U,:).*diams.*log(10),{'Color',rgb('purple'),'LineWidth',2},0.5); hold on
h2=shadedErrorBar(diams,ALL_csa_summer(U,:).*diams.*log(10),ALL_csa_summer_std(U,:).*diams.*log(10),{'Color',rgb('green'),'LineWidth',2},0.5); 
h3=shadedErrorBar(diams,ALL_csa_fall(U,:).*diams.*log(10),ALL_csa_fall_std(U,:).*diams.*log(10),{'Color',rgb('orange'),'LineWidth',2},0.5); 
h4=shadedErrorBar(diams,ALL_csa_winter(U,:).*diams.*log(10),ALL_csa_winter_std(U,:).*diams.*log(10),{'Color',rgb('black'),'LineWidth',2},0.5); 
set(gca,'XScale','log'); axis tight; xlim([3 100])
title('Haptophyta')
xlabel('diamemeter, \mum'); ylabel('\muL L^-^1 per bin')
legend([h1(1).mainLine h2(1).mainLine h3(1).mainLine h4(1).mainLine],'spring','summer','fall','winter')
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 addlabel('f)')
 set(gcf,'Position',[404         185        1225         751])
savefig(gcf,[figsdir filesep 'Fig4.fig']); 
print(gcf,[figsdir filesep 'Fig4'],'-dpng','-r300');




%% Fig 5: Example of how volumetric size distribution of four major taxonomic groups are partitioned among sub-categories during summer
figure(444)
set(gcf,'Position',[68    50   926   946])
subplot(4,3,1)
U = 2; u=[22:27];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'diatom-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 legend('Location','EastOutside')
 ylim([0 0.028])
  vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
  title('Bacillariophyta')
ylabel('\muL L^-^1 per bin')
h= legend('Location','EastOutside');
h.String = {'total','Centric (e)','Hemiaulales (d)','Mastogloia (c)','Pennate (b)','Rhyzosoleniales (a)','other'};

  
  subplot(4,3,4)
U = 9; u=[29:32];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'dinophyceae-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 ylim([0 0.03])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
  title('Dinoflagellata')
 legend('Location','EastOutside')
ylabel('\muL L^-^1 per bin')
h= legend('Location','EastOutside');
h.String = {'total','large (c)','medium (b)','small (a)','other'};

subplot(4,3,7)
U = 6; u=[9:14];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'cyanobacteria-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 ylim([0 0.015])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
  title('Cyanobacteria')
ylabel('\muL L^-^1 per bin'); 
h= legend('Location','EastOutside');
h.String = {'total','clusters (f)','Crocosphaera-like (e)','other (d)','Richelia-free (c)','Trichodesmium-dead (b)','Trichodesmium (a)'};


subplot(4,3,10)
U = 11; u=[34:43];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'haptophyta-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([3 30]);
 ylim([0 0.0075])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
  title('Haptophyta')
 legend('Location','EastOutside')
 % fix name
h= legend('Location','EastOutside');
h.String = {'total','Chrysochromulina-like (i)','Emiliania-like (h)','Discosphaera (g)','Primnesiophyceae-parts (f)','Primnesiophyceae-spherical (e)','Rhabdosphaera (d)','Syracosphaera (c)','Umbellosphaera (b)','Umbilicosphaera (a)','other'};

ylabel('\muL L^-^1 per bin'); xlabel('Diameter [\mum]')

savefig(gcf,[figsdir filesep 'Fig5.fig']); print(gcf,[figsdir filesep 'Fig5'],'-dpng','-r300');


%% Supporting Information Fig. S1: : cruise-averaged volume distributions per month 
load([datadir filesep 'avg_bulk_psd_aloha_kahe_notCSA.mat'],'psdn_aloha','vdn_aloha','gmt_aloha','cruisen_aloha')

u=unique(cruisen_aloha);
ind = find(u==308);
u(ind)=[];
clear XX t
for i = 1:length(u)
    ind = find(cruisen_aloha==u(i))
    XX(i,:)=nanmean(vdn_aloha(ind,:),1); %X(X==0)=NaN;
    yy(i,:) = nanmean(psdn_aloha(ind,:),1);
    t(i,:) = datevec(nanmin(gmt_aloha(ind))); % min or mean will give different results (some cruises started in summer, ended in fall)
end


% one fig per month to make it explicitly clear.
tit = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% force august for t17 (whterver, i did that to force t17 to be summer, but
% ok to show as september here?)
%t(17,2)=8;
figure
for i = 1:12
    subplot(4,3,i)
    ind = find(t(:,2)==i)
    if i==12|i==1|i==2
        colors = 'k';
        leglab = 'Dec, Jan, Feb';
    end
    if i==3|i==4|i==5
        colors=rgb('purple');
        leglab = 'Mar, Apr, May';
    end
      if i==6|i==7|i==8
        colors='g';
        leglab = 'Jun, Jul, Aug';
      end
      if i==9|i==10|i==11
        colors=rgb('orange');
        leglab = 'Sep, Oct, Nov';
      end

      if i==1|i==4|i==7|i==10
          subplot(4,3,i)
        ylabel('\muL L^-^1 per bin')
      end

      if i==10|i==11|i==12
        subplot(4,3,i)
        xlabel('Diameter, \mum')
      end

    semilogx(diams,XX(ind,:).*diams.*log(10),'Color',colors)
    title(tit{i})
    hold on
    ylim([0 0.21]); xlim([3 100]); set(gca,'FontSize',16)
end

subplot(4,3,1)
        ylabel('\muL L^-^1 per bin')
        subplot(4,3,4)
        ylabel('\muL L^-^1 per bin')
        subplot(4,3,7)
        ylabel('\muL L^-^1 per bin')
        subplot(4,3,10)
        ylabel('\muL L^-^1 per bin')
subplot(4,3,10)
   xlabel('Diameter, \mum')
subplot(4,3,11)
   xlabel('Diameter, \mum')
subplot(4,3,12)
   xlabel('Diameter, \mum')

savefig(gcf,[figsdir filesep 'FigS1.fig']); print(gcf,[figsdir filesep 'FigS1'],'-dpng','-r300');


%% Fig. S2: Average seasonal abundances
fig1 = [2 3 5 6 7 9 11 14];
fig2 = [1 4 8 10 12];
fig3 = [13];

class_labels = Labels;

figure
subplot 132
h1=errorbar([1:length(fig1)]-0.1,[abund_spring(fig1)],[abund_spring_se(fig1)],'.','CapSize',0,'Markersize',20); hold on
axis tight; 
h1(1).Color = rgb('purple');
h1=errorbar([1:length(fig1)],[abund_summer(fig1)],[abund_summer_se(fig1)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('green');
h1=errorbar([1:length(fig1)]+0.1,[abund_fall(fig1)],[abund_fall_se(fig1)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('orange');
h1=errorbar([1:length(fig1)]+0.2,[abund_winter(fig1)],[abund_winter_se(fig1)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('black');
set(gca,'XTick',[1:length(fig1)],'XTickLabel',class_labels(fig1))
xtickangle(45)
xlim([-0.1 9.1])

subplot 133
h1=errorbar([1:length(fig2)]-0.1,[abund_spring(fig2)],[abund_spring_se(fig2)],'.','CapSize',0,'Markersize',20); hold on
axis tight; 
h1(1).Color = rgb('purple');
h1=errorbar([1:length(fig2)],[abund_summer(fig2)],[abund_summer_se(fig2)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('green');
h1=errorbar([1:length(fig2)]+0.1,[abund_fall(fig2)],[abund_fall_se(fig2)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('orange');
h1=errorbar([1:length(fig2)]+0.2,[abund_winter(fig2)],[abund_winter_se(fig2)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('black');
set(gca,'XTick',[1:length(fig2)],'XTickLabel',class_labels(fig2))
ylim([0 350])
xtickangle(45)
xlim([0.5 5.6])

subplot 131
h1=errorbar([1:length(fig3)]-0.1,[abund_spring(fig3)],[abund_spring_se(fig3)],'.','CapSize',0,'Markersize',20); hold on
axis tight; 
h1(1).Color = rgb('purple');
h1=errorbar([1:length(fig3)],[abund_summer(fig3)],[abund_summer_se(fig3)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('green');
h1=errorbar([1:length(fig3)]+0.1,[abund_fall(fig3)],[abund_fall_se(fig3)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('orange');
h1=errorbar([1:length(fig3)]+0.2,[abund_winter(fig3)],[abund_winter_se(fig3)],'.','CapSize',0,'Markersize',20)
axis tight; 
h1(1).Color = rgb('black');
set(gca,'XTick',[1:length(fig3)],'XTickLabel',class_labels(fig3))
xlim([0.5 1.5]); ylim([2*10^4 6*10^4 ])
ylabel('particles L^-^1')
xtickangle(45)
legend('spring','summer','fall','winter')
savefig(gcf,[figsdir filesep 'FigS2.fig']); 
print(gcf,[figsdir filesep 'FigS2'],'-dpng','-r300'); % need to modify
%boxes before saving png


%% Fig S4: boxplots of sizes per class 
clear all; 
files = dir('\\jett\awlab\DATA\IFCB CRUISE DATA\HOT\HOT*\proc\CNN_TS4_2025-02-26_Group2_perParticle\');

tax = readtable('\\jett\AWlab\DATA\Fernanda\IFCB\taxonomy_CNN_TS4_2025-02-26.xlsx');
u = unique(tax.GroupsPaper);
c = [1:length(u)];
groups_paper_list = nan(length(c),1);


% assign number to each row in the excel spreadsheet
for i = 1:length(c)
    ind = strmatch(u{i},tax.GroupsPaper);
    groups_paper_list(ind,1)=c(i);
end

% assign a number to each of the files to make it easier to group them
IDX = nan(length(files),1);
for i = 1:length(tax.GroupsPaper)

    ind = strmatch(string(tax.CNN(i)),string({files.name}));

    IDX(ind,1) = groups_paper_list(i);
end

% order ot T_all is same as order of u
clear T_all data
for i = 1:length(c)
    ind = find(IDX==i); % 2 is bad; so will have to get rid of it at some point below
    % load each of the things that  Iwan to merge together, put all into
    % one row
    T=[];
    for j = 1:length(ind)
        load(fullfile(files(ind(j)).folder,files(ind(j)).name));
        data = [ifcbTime lat lon cruise ifcbData.ESD_summed ifcbData.MajorAxisLength ifcbData.MinorAxisLength];
        % choice: remove all but station aloha??? for paper, i think so.
        xxx = find(lat>=22.6 & lat<=22.9 & lon>=-158.2 & lon <=-157.9);
        data=data(xxx,:);
        if ~isempty(data)
            T = [T; data];
        end
    end
    T_all{i,:} = T;
    disp(i)
end




TT=[]; nsamples=[];
for i = 1:length(T_all)
    TT = [TT; T_all{i,1}(:,5) ones(length(T_all{i,1}(:,1)),1)*i];
end

ind = find(TT(:,1)==0) % zeros??? half of "bad" are zeros in diameter. are these skipped in other processing when i look for >3-98?
TT(ind,:)=[];
ind = find(TT(:,2)==2) % bad
TT(ind,:)=[];

clear nsamples
for i = 1:length(c)
    ind= find(TT(:,2)==i);
    nsamples(i,:) = length(ind);
end

nsamples(2)=[];
u(2)=[];
c(2)=[];

% since i removed "bad", which was colum=2; i need to change all numbers (otherwise there's a gap in boxplots between artropoda and chloro but labels are mismatched):
ind=find(TT(:,2)>1);
TT(ind,2)=TT(ind,2)-1;


% change the labels for some:
u{19-1,1}='dictyochophyta';
u{20-1,1}='dinoflagellata-large';
u{21-1,1}='dinoflagellata-medium';
u{22-1,1}='dinoflagellata-small';
u{23-1,1}='dinoflagellata-unknown';
u{25,1}='haptophyta-emiliania-like';
u{24,1}='haptophyta-chrysochromulina-like';
u{34-1,1}='haptophyta-unknown';

% append number of samples 
labels=strcat(string(u),' (',num2str(nsamples),')')

 load('\\jett\awlab\DATA\Fernanda\PAPERS\HOT_PSD\submission\code\data\max_good_diam.mat')

figure
boxchart(TT(:,2),TT(:,1),'Orientation','horizontal','MarkerStyle','none','MarkerSize',1); hold on
set(gca,'YTick',[1:36],'YTickLabel',labels)
set(gca,'YDir','reverse')
ylim([0 37])
xlim([0 100])
vline(nanmean(max_good_diam),'--k')
vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
xlabel('Diameter (ESD), \mum')

savefig(gcf,[figsdir filesep 'FigS4.fig']); print(gcf,[figsdir filesep 'FigS4'],'-dpng','-r300');

% major and minor
TT=[]; nsamples=[];
for i = 1:length(T_all)
    TT = [TT; T_all{i,1}(:,6) T_all{i,1}(:,7) ones(length(T_all{i,1}(:,1)),1)*i];
end

ind = find(TT(:,1)==0) % zeros??? half of "bad" are zeros in diameter. are these skipped in other processing when i look for >3-98?
TT(ind,:)=[];
ind = find(TT(:,2)==2) % bad
TT(ind,:)=[];

clear nsamples
for i = 1:length(c)
    ind= find(TT(:,2)==i);
    nsamples(i,:) = length(ind);
end

nsamples(2)=[];
u(2)=[];
c(2)=[];

% since i removed "bad", which was colum=2; i need to change all numbers (otherwise there's a gap in boxplots between artropoda and chloro but labels are mismatched):
ind=find(TT(:,2)>1);
TT(ind,2)=TT(ind,2)-1;

figure
subplot 121
boxchart(TT(:,3),TT(:,2),'Orientation','horizontal','MarkerStyle','none','MarkerSize',1); hold on
set(gca,'YTick',[1:36],'YTickLabel',labels)
set(gca,'YDir','reverse')
ylim([0 37])
xlim([0 250])
vline(nanmean(max_good_diam),'--k')
vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
xlabel('Minor Axis Length, \mum')

subplot 122
boxchart(TT(:,3),TT(:,1),'Orientation','horizontal','MarkerStyle','none','MarkerSize',1); hold on
set(gca,'YTick',[1:36],'YTickLabel',labels)
set(gca,'YDir','reverse')
ylim([0 37])
xlim([0 250])
vline(nanmean(max_good_diam),'--k')
vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
xlabel('Major Axis Length, \mum')



%% Supporting Information S5 (same as Fig 6 but all seasons)
figure(4441)
set(gcf,'Position',[68    50   926   946])

subplot(4,4,1)
U = 2; u=[22:27];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'diatom-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 legend('Location','EastOutside')
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
  vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
ylabel('\muL L^-^1 per bin')

  
  subplot(4,4,5)
U = 9; u=[29:32];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'dinophyceae-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 legend('Location','EastOutside')
ylabel('\muL L^-^1 per bin')

subplot(4,4,9)
U = 6; u=[9:14];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'cyanobacteria-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 legend('Location','EastOutside')
ylabel('\muL L^-^1 per bin'); 

subplot(4,4,13)
U = 11; u=[34:43];
area(diams,ALL_csa_summer(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_summer_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'haptophyta-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([3 30]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 legend('Location','EastOutside')
ylabel('\muL L^-^1 per bin'); xlabel('Diameter [\mum]')


% add fall
subplot(4,4,2)
U = 2; u=[22:27];
area(diams,ALL_csa_fall(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_fall_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'diatom-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
  vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

  
  subplot(4,4,6)
U = 9; u=[29:32];
area(diams,ALL_csa_fall(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_fall_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'dinophyceae-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

subplot(4,4,10)
U = 6; u=[9:14];
area(diams,ALL_csa_fall(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_fall_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'cyanobacteria-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

subplot(4,4,14)
U = 11; u=[34:43];
area(diams,ALL_csa_fall(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_fall_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'haptophyta-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([3 30]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 xlabel('Diameter [\mum]')


% add winter
subplot(4,4,3)
U = 2; u=[22:27];
area(diams,ALL_csa_winter(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_winter_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'diatom-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
  vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

  
  subplot(4,4,7)
U = 9; u=[29:32];
area(diams,ALL_csa_winter(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_winter_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'dinophyceae-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

subplot(4,4,11)
U = 6; u=[9:14];
area(diams,ALL_csa_winter(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_winter_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'cyanobacteria-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

subplot(4,4,15)
U = 11; u=[34:43];
area(diams,ALL_csa_winter(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_winter_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'haptophyta-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([3 30]);
  xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
 xlabel('Diameter [\mum]')

% add spring
subplot(4,4,4)
U = 2; u=[22:27];
area(diams,ALL_csa_spring(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_spring_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'diatom-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
  vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')

  
  subplot(4,4,8)
U = 9; u=[29:32];
area(diams,ALL_csa_spring(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_spring_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'dinophyceae-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')


subplot(4,4,12)
U = 6; u=[9:14];
area(diams,ALL_csa_spring(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_spring_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'cyanobacteria-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([5 100]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')


subplot(4,4,16)
U = 11; u=[34:43];
area(diams,ALL_csa_spring(U,:).*diams.*log(10),'FaceColor',[.7 .7 .7],'DisplayName','total');
hold on
set (gca, 'Xscale', 'log')
colors = cmocean('balance',length(u)); % to make it match original colors when radiozoa et cwere not removed
for i = 1:length(u)
    plot(diams,ALL_csa_spring_g2(u(i),:).*diams.*log(10),'Color',colors(i,:),'LineWidth',3,'DisplayName',string(extractBetween(f2(u(i)).name,'haptophyta-','_uw_all.mat'))); hold on   
end
 axis tight; xlim([3 30]);
 xx = nanmax([ALL_csa_summer(U,:).*diams.*log(10) ALL_csa_spring(U,:).*diams.*log(10)  ALL_csa_fall(U,:).*diams.*log(10) ALL_csa_winter(U,:).*diams.*log(10) ])
 ylim([0 xx])
 vline(nanmean(max_good_diam),'--k')
 vline(nanmean(max_good_diam)-nanstd(max_good_diam),':k')
 vline(nanmean(max_good_diam)+nanstd(max_good_diam),':k')
xlabel('Diameter [\mum]')

savefig(gcf,[figsdir filesep 'FigS5.fig']); print(gcf,[figsdir filesep 'FigS5'],'-dpng','-r300');




%% RUN STATS (for Supporting Information Tables S2-S5)
% evaluate seasonality (Supporting Information Table S2)
FORPAPER_anova_test_seasonality
% regressions against SD50 and D90 (Supporting Information Table S3 and S4)
runTaxonMomentsStats_MAIN
% structural shape sensitivity test (Supporting Information Table S5)
FORPAPER_mean_shape_contribution_tests

%% get other paper stats

% percent volume contained above and below 20 um pe season
% proportion of below 20um
sum(trapz(diams(1:27),ALL_csa_spring(:,[1:27]),2))./sum(trapz(diams,ALL_csa_spring,2)) % 70% --> 79%
sum(trapz(diams(1:27),ALL_csa_summer(:,[1:27]),2))./sum(trapz(diams,ALL_csa_summer,2)) % 41% --> 60
sum(trapz(diams(1:27),ALL_csa_winter(:,[1:27]),2))./sum(trapz(diams,ALL_csa_winter,2)) % 67% --> 76
sum(trapz(diams(1:27),ALL_csa_fall(:,[1:27]),2))./sum(trapz(diams,ALL_csa_fall,2)) % 54% --> 70

% percent of total per category:
tot = nansum(data,2);
pct = data./tot*100;

% percent particle volume dominated by detritus above and below 20um

sum(trapz(diams(28:end),ALL_csa_summer(7,[28:end]),2))./sum(trapz(diams,ALL_csa_summer(7,:),2)) % 34%
sum(trapz(diams(28:end),ALL_csa_winter(7,[28:end]),2))./sum(trapz(diams,ALL_csa_winter(7,:),2)) % 37%
sum(trapz(diams(28:end),ALL_csa_spring(7,[28:end]),2))./sum(trapz(diams,ALL_csa_spring(7,:),2)) % 57%
sum(trapz(diams(28:end),ALL_csa_fall(7,[28:end]),2))./sum(trapz(diams,ALL_csa_fall(7,:),2)) % 36%

% correlation between diatoms and total volume
for i = 1:14;
    r=corrcoef(tot_volconc,ALL_volconc(i,:)) ;
    rr(i,:)=r(1,2).^2;
end
figure
plot([1:14],rr); set(gca,'XTick',[1:14],'XTickLabel',Labels); % high r with diatoms, dinos, unknown

% diatoms vs cyano:
[r,p]=corrcoef(ALL_volconc(2,:),ALL_volconc(6,:)); % r=0.17, p=0.24

% diatoms vs detritus:
[r,p]=corrcoef(ALL_volconc(2,:),ALL_volconc(7,:)) % r=0.13, p=0.24



% what percent of cyano is tricho and trichdead
trapz(diams,ALL_csa_fall_g2(14,:),2)./ trapz(diams,ALL_csa_fall(6,:),2)% 85
trapz(diams,ALL_csa_fall_g2(13,:),2)./ trapz(diams,ALL_csa_fall(6,:),2)% +2, so 87% is tricho. 




%% Make csv files for Zenodo: major category psds (volume) per cruise

clearvars -except max_good_diam*  good_cruises figsdir nsamples_percruise datadir
f = dir([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1'])
f([1:2 13 15])=[];
%f([1:2])=[];
Labels = extractBefore({f.name},'_uw_all.mat'); Labels{1,8}='Dictyochophyta'; Labels{1,9}='Dinoflagellata'; Labels{1,2}='Bacillariophyta';

good_cruises = [1:4 6:45]; % remove 308

load([datadir filesep 'binWidth']);

Tall_phylum=[];
for i = 1:length(f)
    load(fullfile(f(i).folder,f(i).name),'vdn_aloha','diams','gmt_all','cruisen_all');
    
    vdn_aloha=round(vdn_aloha(good_cruises,:),7);
    gmt_all=gmt_all(good_cruises);
    cruise = cruisen_all(good_cruises);

    
    A = [gmt_all cruise vdn_aloha];
        % convert to table
    T = array2table(A);
    T.A1=datetime(T.A1,'ConvertFrom','datenum');
    
    % add group column
    T.group = repmat(string(Labels{i}), height(T), 1);

    % rearrange
    T = movevars(T,'group','After','A2');

    % append
    Tall_phylum = [Tall_phylum; T];

end

% write csv
writetable(Tall_phylum,'major_categories_PSD_volume_IFCB_HOT294_to_HOT361.csv')


outfile = 'major_categories_PSD_volume_IFCB_HOT294_to_HOT361.csv';

fid = fopen(outfile,'w');

str = strjoin(compose('%.2f', diams), ', ');
str2 = strjoin(compose('%.2f', binWidth), ', ');
str3 = strjoin(compose('HOT%d (%d)', cruise, nsamples_percruise), '; ');

% write custom header lines
fprintf(fid,'Volumetric particle size distributions for 14 major IFCB categories which make up the bulk IFCB particle volume at Station ALOHA averaged per cruise.\n');
fprintf(fid,['Number of IFCB samples (each sample quantifies particles in the ~3-100 micron size range contained  in ~4mL of seawater) averaged per cruise:' str3 '\n']);
fprintf(fid,'Columns: time (average cruise date/time, UTC), HOT cruise number, annotation group (from Convolution Neural Network), columns 4-53 are volume concentration (microliters per liter per micron, with volume calculated from biovolume-derived equivalent spherical diameter) for 50 log-spaced diameter bins.\n');
fprintf(fid,['Diameter bin centers (micrometers):' str '\n']);
fprintf(fid,['Diameter binwidths (micrometers):' str2 '\n']);
fprintf(fid,'Note: integrating volume concentration over the 50 diameter bins yields total volume concentration in units of microliters per liter per time step and taxa. Summing all volume concentrations per taxa for a single time step yields total volume concentration imaged by the IFCB (in the 3-100 micron size range). \n');
fprintf(fid,'F1-scores quantifying accuracy of particle annotations for the different groups are provided separately.\n');
fprintf(fid,'IFCB Dashboard is available here: https://ifcbdb2.soest.hawaii.edu/datasets\n');
fprintf(fid,'Generated on: %s\n', datestr(now));
fprintf(fid,'\n');   % blank line

fclose(fid);

% append the table after the header
writetable(Tall_phylum, outfile, 'WriteMode','append');


%% Make csv files for Zenodo: major category psds (abundances) per cruise
clearvars -except max_good_diam*  good_cruises figsdir nsamples_percruise datadir
f = dir([datadir filesep 'CNN_TS4_2025-02-26_merged_Group1'])
f([1:2 13 15])=[];
%f([1:2])=[];
Labels = extractBefore({f.name},'_uw_all.mat'); Labels{1,8}='Dictyochophyta'; Labels{1,9}='Dinoflagellata'; Labels{1,2}='Bacillariophyta';

good_cruises = [1:4 6:45]; % remove 308

load([datadir filesep 'binWidth']);

Tall_phylum=[];
for i = 1:length(f)
    load(fullfile(f(i).folder,f(i).name),'psdn_aloha','diams','gmt_all','cruisen_all');
    
    p_aloha=round(psdn_aloha(good_cruises,:),4);
    gmt_all=gmt_all(good_cruises);
    cruise = cruisen_all(good_cruises);

    
    A = [gmt_all cruise p_aloha];
        % convert to table
    T = array2table(A);
    T.A1=datetime(T.A1,'ConvertFrom','datenum');
    
    % add group column
    T.group = repmat(string(Labels{i}), height(T), 1);

    % rearrange
    T = movevars(T,'group','After','A2');

    % append
    Tall_phylum = [Tall_phylum; T];

end

% write csv
writetable(Tall_phylum,'major_categories_PSD_abundance_IFCB_HOT294_to_HOT361.csv')


outfile = 'major_categories_PSD_abundance_IFCB_HOT294_to_HOT361.csv';

fid = fopen(outfile,'w');

str = strjoin(compose('%.2f', diams), ', ');
str2 = strjoin(compose('%.2f', binWidth), ', ');
str3 = strjoin(compose('HOT%d (%d)', cruise, nsamples_percruise), '; ');

% write custom header lines
fprintf(fid,'Particle size distributions (particles per liter per micron) for 14 major IFCB categories which make up the bulk IFCB total particle numbers at Station ALOHA averaged per cruise.\n');
fprintf(fid,['Number of IFCB samples (each sample quantifies particles in the ~3-100 micron size range contained  in ~4mL of seawater) averaged per cruise:' str3 '\n']);
fprintf(fid,'Columns: time (average cruise date/time, UTC), HOT cruise number, annotation group (from Convolution Neural Network), columns 4-53 are particle abundances or counts (number of particles per liter per micron) for 50 log-spaced diameter bins.\n');
fprintf(fid,['Diameter bin centers (micrometers):' str '\n']);
fprintf(fid,['Diameter binwidths (micrometers):' str2 '\n']);
fprintf(fid,'Note: integrating particle counts over the 50 diameter bins yields total number of particles in per liter units per time step and taxa. Summing all counts per liter per taxa for a single time step yields total number of images imaged by the IFCB (in the 3-100 micron size range). \n');
fprintf(fid,'F1-scores quantifying accuracy of particle annotations for the different groups are provided separately.\n');
fprintf(fid,'IFCB Dashboard is available here: https://ifcbdb2.soest.hawaii.edu/datasets\n');
fprintf(fid,'Generated on: %s\n', datestr(now));
fprintf(fid,'\n');   % blank line

fclose(fid);

% append the table after the header
writetable(Tall_phylum, outfile, 'WriteMode','append');



