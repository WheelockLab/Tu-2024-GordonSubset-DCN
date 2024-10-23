clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
adult_color = [236,0,138]/255; % pink, Figure 1A
infant_color = [127,63,152]/255; % purple, Figure 1B

clim = [-0.5,0.5]; % color scale for FC

T = readtable('BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_stats_withMullen.csv');

load('BCP_Gordon_BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_7pt2min_randsample.mat')
zmatBCP = zmat;
for ii = 1:size(zmatBCP,1),zmatBCP(ii,ii,:) = 0;end
avg_zmatBCP = mean(zmatBCP,3);

load('/data/wheelock/data1/datasets/WashU120/pconns/washu120_parcellation_Gordon_20231101.mat')
for ii = 1:size(zmat,1),  zmat(ii,ii,:) = 0; end
zmatWashU120 = mean(zmat,3);

load('IM_Gordon_13nets_333Parcels_renamed.mat','IM')
zmat_Gordon_all120 = zmat(IM.order,IM.order,:);
zmat_gordon_BCP = zmatBCP(IM.order,IM.order,:);
zmat_gordon_WashU120 = zmatWashU120(IM.order,IM.order);

avg_zmat_gordon_BCP = avg_zmatBCP(IM.order,IM.order);

KardanIM = load('IM_11_BCP94_renamed.mat','IM');
KardanIM = KardanIM.IM;
zmat_Kardan_BCP = zmatBCP(KardanIM.order,KardanIM.order,:);
zmat_Kardan_WashU120 = zmatWashU120(KardanIM.order,KardanIM.order);
zmat_Kardan_all120 = zmat(KardanIM.order,KardanIM.order,:);

avg_zmat_Kardan_BCP = avg_zmatBCP(KardanIM.order,KardanIM.order);%mean(zmat_Kardan_BCP,3);

[~,Nroi,~]=size(zmat_gordon_BCP);
%% Calculate SI for Adult and Infant FC
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));

D = calc_correlationdist(zmat_gordon_WashU120);
s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
SI_Adult_Gordon = NaN(Nroi,1);
SI_Adult_Gordon(keepnets) = s;

D = calc_correlationdist(avg_zmat_gordon_BCP);
s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
SI_Infant_Gordon = NaN(Nroi,1);
SI_Infant_Gordon(keepnets) = s;

noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
keepnets = KardanIM.key(:,2)~=noneidx;
M = ones(max(KardanIM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));

D = calc_correlationdist(avg_zmat_Kardan_BCP);
s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
SI_Infant_Kardan = NaN(Nroi,1);
SI_Infant_Kardan(keepnets) = s;

D = calc_correlationdist(zmat_Kardan_WashU120);
s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
SI_Adult_Kardan = NaN(Nroi,1);
SI_Adult_Kardan(keepnets) = s;
%% A big confound is that whether we get high SI across the group because they are the locations with low interindividual variability in infants, 
% so we calculate SI on individuals here
Nsess = size(zmat_Gordon_all120,3);
SI_YA_Gordon_session = NaN(Nroi,Nsess);
SI_YA_Kardan_session = NaN(Nroi,Nsess);

for isess = 1:Nsess
    isess
    noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
    keepnets = IM.key(:,2)~=noneidx;
    M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    D = calc_correlationdist(zmat_Gordon_all120(:,:,isess));
    s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
    SI_YA_Gordon_session(keepnets,isess) = s;
    
    noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
    keepnets = KardanIM.key(:,2)~=noneidx;
    M = ones(max(KardanIM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    D = calc_correlationdist(zmat_Kardan_all120(:,:,isess));
    s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
    SI_YA_Kardan_session(keepnets,isess) = s;
end

%% Supplementary: visualize the SI across individual sessions (WU120)
cmap = interp1(linspace(0,100,10),redbluecmap(10),linspace(0,100,100)); % red blue cmap but interploated
sortid = 1:size(SI_YA_Gordon_session,2);

figure('units','inches','position',[1,1,3,3]);
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
netcmap = [IM.cMap;1,1,1];
key = {IM.key(keepnets,:),repelem(14,Nsess,2)};
Matrix_Org3(SI_YA_Gordon_session(keepnets,sortid),key,10,[-1,1],netcmap,0,cmap)
xlabel('Sessions');
ylabel('Areas');
% title('Adult Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Gordon.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);
noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
keepnets = KardanIM.key(:,2)~=noneidx;
netcmap = [KardanIM.cMap;1,1,1];
key = {KardanIM.key(keepnets,:),repelem(12,Nsess,2)};
Matrix_Org3(SI_YA_Kardan_session(keepnets,sortid),key,10,[-1,1],netcmap,0,cmap)
xlabel('Sessions');
ylabel('Areas');
% title('Infant Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Kardan.tif','-dtiff','-r300');

makecolorbar(cmap,[-1,1],'h','SI')
% print('./Figures/SI_color.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,2,3]);hold on
X = [nanmean(SI_YA_Gordon_session);nanmean(SI_YA_Kardan_session)]';
for isess = 1:Nsess
    plot([0,1],X(isess,:),'Color',[0.5,0.5,0.5,0.1]);
end
boxplot(X,1:2,'position',[0,1],'colors',[adult_color;infant_color],'symbol','.'); 
h = findobj(gca, 'Tag', 'Box');h2 = findobj(gca,'Tag','Upper Whisker');h3 = findobj(gca,'Tag','Lower Whisker');
for k =  1:length(h)
    h(k).LineWidth = 2; h2(k).LineWidth = 2; h3(k).LineWidth = 2;
end
[h,p] = ttest(X(:,1),X(:,2))
sigstar({[0,1]},p);
xticklabels({'Adult Networks','Infant Networks'})
xtickangle(45);
ylabel('avg SI ');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_pairedttest.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,1.5,1.5]);hold on;
[r,p] = corr(X(:,1),X(:,2));
scatter(X(:,1),X(:,2),15,'k','.')
xl = xlim; yl = ylim;
lims = [min([xl(1),yl(1)]),max([xl(2),yl(2)])];
xlim(lims);ylim(lims);
axis equal;axis square;
hl =refline;hl.Color= 'r';
xlim(lims);ylim(lims);
text(0.05,0.8,sprintf('r = %1.2f',r),'units','normalized')
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_GordonKardanScatter.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);hold on;
xlim([-1,1]);ylim([-1,1]);
xline(0,'LineStyle','-');yline(0,'LineStyle','-');
errorbar(SI_Adult_Gordon,mean(SI_YA_Gordon_session,2),std(SI_YA_Gordon_session,[],2),'LineStyle','None','CapSize',1,'Color',[0.8,0.8,0.8,0.1]);
scatter(SI_Adult_Gordon,mean(SI_YA_Gordon_session,2),15,adult_color,'.');
h = refline;h.Color = 'r';
[r,p] = corr(SI_Adult_Gordon,mean(SI_YA_Gordon_session,2),'rows','complete')
text(0.65,0.1,sprintf('r = %1.2f',r),'Units','normalized','FontSize',12);
xlabel('SI on group average FC');
ylabel('SI on individual FC');
title('Adult Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Gordon_GroupVSIndividual.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);hold on;
xlim([-1,1]);ylim([-1,1]);
xline(0,'LineStyle','-');yline(0,'LineStyle','-');
errorbar(SI_Adult_Kardan,mean(SI_YA_Kardan_session,2),std(SI_YA_Kardan_session,[],2),'LineStyle','None','CapSize',1,'Color',[0.8,0.8,0.8,0.1]);
scatter(SI_Adult_Kardan,mean(SI_YA_Kardan_session,2),15,infant_color,'.');
h = refline;h.Color = 'r';
[r,p] = corr(SI_Adult_Kardan,mean(SI_YA_Kardan_session,2),'rows','complete')
text(0.65,0.1,sprintf('r = %1.2f',r),'Units','normalized','FontSize',12);
xlabel('SI on group average FC');
ylabel('SI on individual FC');
title('Infant Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Kardan_GroupVSIndividual.tif','-dtiff','-r300');

% some stats
[h,p] = ttest(X,0) % compare individual session SI to 0

mean(X)
std(X)

mean(X(:,1)-X(:,2))
std(X(:,1)-X(:,2))
computeCohen_d(X(:,1),X(:,2),'paired')

h = arrayfun(@(iROI)ttest(SI_YA_Gordon_session(iROI,:),0,'tail','right'),1:333); nansum(h)
h = arrayfun(@(iROI)ttest(SI_YA_Kardan_session(iROI,:),0,'tail','right'),1:333); nansum(h)
%% Supplementary: visualize the SI across individual sessions (BCP)
Nsess = size(zmat_gordon_BCP,3);
SI_BCP_Gordon_session = NaN(333,Nsess);
SI_BCP_Kardan_session = NaN(333,Nsess);
for isess = 1:Nsess
    isess
    noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
    keepnets = IM.key(:,2)~=noneidx;
    M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    D = calc_correlationdist(zmat_gordon_BCP(:,:,isess));
    s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
    SI_BCP_Gordon_session(keepnets,isess) = s;
    
    noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
    keepnets = KardanIM.key(:,2)~=noneidx;
    M = ones(max(KardanIM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    D = calc_correlationdist(zmat_Kardan_BCP(:,:,isess));
    s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
    SI_BCP_Kardan_session(keepnets,isess) = s;
end

cmap = interp1(linspace(0,100,10),redbluecmap(10),linspace(0,100,100)); % red blue cmap but interploated
[~,sortid] = sort(T.age_yrs);
%% Plot (BCP)
figure('units','inches','position',[1,1,3,3]);
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
netcmap = [IM.cMap;1,1,1];
key = {IM.key(keepnets,:),repelem(14,Nsess,2)};
Matrix_Org3(SI_BCP_Gordon_session(keepnets,sortid),key,10,[-1,1],netcmap,0,cmap)
xlabel('Sessions');
ylabel('Areas');
title('Adult Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Gordon.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);
noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
keepnets = KardanIM.key(:,2)~=noneidx;
netcmap = [KardanIM.cMap;1,1,1];
key = {KardanIM.key(keepnets,:),repelem(12,Nsess,2)};
Matrix_Org3(SI_BCP_Kardan_session(keepnets,sortid),key,10,[-1,1],netcmap,0,cmap)
xlabel('Sessions');
ylabel('Areas');
title('Infant Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Kardan.tif','-dtiff','-r300');

makecolorbar(cmap,[-1,1],'h','SI')
% print('./Figures/SI_color.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,2,3]);hold on
X = [nanmean(SI_BCP_Gordon_session);nanmean(SI_BCP_Kardan_session)]';
for isess = 1:Nsess
    plot([0,1],X(isess,:),'Color',[0.5,0.5,0.5,0.1]);
end
boxplot(X,1:2,'position',[0,1],'colors',[adult_color;infant_color],'symbol','.'); 
h = findobj(gca, 'Tag', 'Box');h2 = findobj(gca,'Tag','Upper Whisker');h3 = findobj(gca,'Tag','Lower Whisker');
for k =  1:length(h)
    h(k).LineWidth = 2; h2(k).LineWidth = 2; h3(k).LineWidth = 2;
end
[h,p] = ttest(X(:,1),X(:,2))
sigstar({[0,1]},p);
xticklabels({'Adult Networks','Infant Networks'});
xtickangle(45)
ylabel('avg SI ');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_pairedttest.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,1.5,1.5]);hold on;
[r,p] = corr(X(:,1),X(:,2));
scatter(X(:,1),X(:,2),15,'k','.')
xl = xlim; yl = ylim;
lims = [min([xl(1),yl(1)]),max([xl(2),yl(2)])];
xlim(lims);ylim(lims);
axis equal;axis square;
hl =refline;hl.Color= 'r';
xlim(lims);ylim(lims);
text(0.5,0.05,sprintf('r = %1.2f',r),'units','normalized')
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_GordonKardanScatter.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);hold on;
xlim([-1,1]);ylim([-1,1]);
xline(0,'LineStyle','-');yline(0,'LineStyle','-');
errorbar(SI_Infant_Gordon,mean(SI_BCP_Gordon_session,2),std(SI_BCP_Gordon_session,[],2),'LineStyle','None','CapSize',1,'Color',[0.8,0.8,0.8,0.1]);
scatter(SI_Infant_Gordon,mean(SI_BCP_Gordon_session,2),15,adult_color,'.');
h = refline;h.Color = 'r';
[r,p] = corr(SI_Infant_Gordon,mean(SI_BCP_Gordon_session,2),'rows','complete')
text(0.65,0.1,sprintf('r = %1.2f',r),'Units','normalized','FontSize',12);
xlabel('SI on group average FC');
ylabel('SI on individual FC');
title('Adult Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Gordon_GroupVSIndividual.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);hold on;
xlim([-1,1]);ylim([-1,1]);
xline(0,'LineStyle','-');yline(0,'LineStyle','-');
errorbar(SI_Infant_Kardan,mean(SI_BCP_Kardan_session,2),std(SI_BCP_Kardan_session,[],2),'LineStyle','None','CapSize',1,'Color',[0.8,0.8,0.8,0.1]);
scatter(SI_Infant_Kardan,mean(SI_BCP_Kardan_session,2),15,infant_color,'.');
h = refline;h.Color = 'r';
[r,p] = corr(SI_Infant_Kardan,mean(SI_BCP_Kardan_session,2),'rows','complete')
text(0.65,0.1,sprintf('r = %1.2f',r),'Units','normalized','FontSize',12);
xlabel('SI on group average FC');
ylabel('SI on individual FC');
title('Infant Networks');
set(gca,'FontSize',12);
% print('./Figures/IndividualSI_Kardan_GroupVSIndividual.tif','-dtiff','-r300');

% some stats
[h,p] = ttest(X,0) % compare individual session SI to 0

mean(X)
std(X)

mean(X(:,1)-X(:,2))
std(X(:,1)-X(:,2))
computeCohen_d(X(:,1),X(:,2),'paired')

h = arrayfun(@(iROI)ttest(SI_BCP_Gordon_session(iROI,:),0,'tail','right'),1:333); nansum(h)
h = arrayfun(@(iROI)ttest(SI_BCP_Kardan_session(iROI,:),0,'tail','right'),1:333); nansum(h)
% Plot colorbar
[hCB,hf] = makecolorbar(jet,clim,'v','z(r)')

%% Is there a correlation between age and SI?

figure('units','inches','position',[1,1,3,3]);hold on;
scatter(T.age_yrs,nanmean(SI_BCP_Kardan_session),15,'k.');
[r,p] = corr(T.age_yrs,[nanmean(SI_BCP_Kardan_session)]')
h = refline; h.Color = 'r';
xlabel('Age (yrs)');
ylabel({'avg SI';'(Kardan)'});
set(gca,'FontSize',12);
print('./Figures/IndividualSI_Kardan_Infant.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);hold on;
scatter(T.age_yrs,nanmean(SI_BCP_Gordon_session),15,'k.');
[r,p] = corr(T.age_yrs,[nanmean(SI_BCP_Gordon_session)]')
h = refline; h.Color = 'r';
xlabel('Age (yrs)');
ylabel({'avg SI';'(Gordon)'});
set(gca,'FontSize',12);
print('./Figures/IndividualSI_Gordon_Infant.tif','-dtiff','-r300');

figure('units','inches','position',[1,1,3,3]);hold on;
scatter(T.age_yrs,nanmean(SI_BCP_Kardan_session)-nanmean(SI_BCP_Gordon_session),15,'k.');
[r,p] = corr(T.age_yrs,[nanmean(SI_BCP_Kardan_session)-nanmean(SI_BCP_Gordon_session)]')
h = refline; h.Color = 'r';
xlabel('Age (yrs)');
ylabel({'avg SI difference';'(Kardan-Gordon)'});
set(gca,'FontSize',12);
print('./Figures/IndividualSI_diff_Infant.tif','-dtiff','-r300');
