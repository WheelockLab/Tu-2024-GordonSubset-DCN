clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
T = readtable('BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_stats_withMullen.csv');

% Load Infant
datafile = 'BCP_Tu_326_BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_7pt2min_randsample.mat'
load(datafile);
zmatBCP = zmat;
for ii = 1:size(zmatBCP,1),zmatBCP(ii,ii,:) = 0;end % make diagonals zero
avg_zmatBCP = mean(zmatBCP,3);

[~,Nroi]=size(avg_zmatBCP);
clim = [-0.5,0.5]; % color scale for FC
load('Parcels_Tu_326.mat');
network_name = '12Networks'
load(['IM_Tu326_',network_name,'.mat'],'IM')
avg_zmatBCP_Tu326 = avg_zmatBCP(IM.order,IM.order);
zmatBCP_Tu326 = zmatBCP(IM.order,IM.order,:);
%% ============== Repeat Figure 1 ====================
%% Plot FC sorted by Tu 12 Networks
keepnets = true(Nroi,1);
M = ones(max(IM.key(:,2)));M = M-diag(diag(M));

figure;
Matrix_Org3(avg_zmatBCP_Tu326(keepnets,keepnets),...
    IM.key(keepnets,:),10,clim,IM.cMap,0);
title('Infant FC');
D = calc_correlationdist(avg_zmatBCP_Tu326);
s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
text(0.2,-0.1,sprintf('avg SI = %2.3f',mean(s)),'Units','normalized','FontWeight','Bold');
print(['./Figures/FcOrganizedByTu326Network'],'-dpng');
%% Plot values on brain
[~,sortid] = sort(IM.order);
SI_Infant_Tu326 = s(sortid);
cmap = interp1(linspace(0,100,10),redbluecmap(10),linspace(0,100,100)); % red blue cmap but interploated

figure;
tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot_parcels_by_values(SI_Infant_Tu326,'med',Parcels,[-1,1],cmap) 
nexttile;
plot_parcels_by_values(SI_Infant_Tu326,'lat',Parcels,[-1,1],cmap) 
print(['./Figures/FcOrganizedByTu326Network_glassbrain'],'-dpng');
%% try to do bootstrap CI for the average SI estimate for BCP
rng('default')
SilTu326Mean = NaN(1000,1);
N = size(zmatBCP_Tu326,3)
for isample = 1:1000
    isample
    idx = randsample(N,N,1);
    
    zmatBCP_Tu326_subsample = mean(zmatBCP_Tu326(:,:,idx),3);
    D = calc_correlationdist(zmatBCP_Tu326_subsample);
    keepidx = true(Nroi,1);
    M = ones(max(IM.key(:,2)));M = M-diag(diag(M));

    s = silhouette_coef_mod(IM.key(keepidx,2),D(keepidx,keepidx),M);
    SilTu326Mean(isample) = mean(s);    
end

CI95 = quantile(SilTu326Mean,[0.025,0.975])
%% ============== Repeat Figure 4 ====================
%% Moving average window for age
noneidx = [];
keepnets = true(Nroi,1);
M = ones(max(IM.key(:,2)));M = M-diag(diag(M));

[Ages,AgesortID] = sort(T.age_yrs);
Ages = Ages(1:281);AgesortID = AgesortID(1:281); % after 281 the data point is very scattered
windowsz = 20;% 20 sessions
windowstep = 1;
N = floor((length(Ages)-windowsz)/windowstep);
[sil_Tu326,agemean] = deal(NaN(N,1));

for jj =1:N
    jj
    agemean(jj) =mean(Ages([1:windowsz]+(jj-1)*windowstep)); % ages weighted by the mean
    wholesample = AgesortID([1:windowsz]+(jj-1)*windowstep);
    tmpzmat=mean(zmatBCP_Tu326(:,:,wholesample),3);
    D =calc_correlationdist(tmpzmat);
    s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
    sil_Tu326(jj) = mean(s);
end
save(['./Results/moving_avg_results_Tu326_',network_name,'.mat'],'sil_Tu326*','agemean');
%%
Tu326_12Networks = load('./Results/moving_avg_results_Tu326_12Networks.mat','sil_Tu326')
Tu326_19Networks = load('./Results/moving_avg_results_Tu326_19Networks.mat','sil_Tu326')
load('moving_avg_results.mat')
legendstr = {'Gordon','Gordon (Subset)','Kardan','Kardan (Subset)','Tu326_12Networks','Tu326_19Networks'}
clear h
clear h
figure('position',[100 100 400 400]);hold on;
h(1) = plot(agemean,sil_Gordon,'LineWidth',2);
h(2) = plot(agemean,sil_Gordon_subset,'LineWidth',2);
h(3) = plot(agemean,sil_Kardan,'LineWidth',2);
h(4) = plot(agemean,sil_Kardan_subset,'LineWidth',2);
h(5) = plot(agemean,Tu326_12Networks.sil_Tu326,'LineWidth',2);
h(6) = plot(agemean,Tu326_19Networks.sil_Tu326,'LineWidth',2);
xlim([0.5,2.5]);
Plot.vline(agemean(53)); %~ 1 year
Plot.vline(agemean(223)); % ~2  year
legend(h,legendstr,'location','best');
legend('location','southoutside');
ylabel('SI');
xlabel('Age (yrs)')
set(gca,'FontSize',12);
xlim([0.5,2.5]);
ylim([-0.1,0.4]);
xline(agemean(53),'LineStyle','--'); %~ 1 year
xline(agemean(223),'LineStyle','--'); % ~2  year
legend(h,legendstr,'location','best','interpreter','None');
legend('location','southoutside');
ylabel('SI');
xlabel('Age (yrs)')
set(gca,'FontSize',12);

print(['./Figures/SilAgeWindowTu326Network'],'-dpng');