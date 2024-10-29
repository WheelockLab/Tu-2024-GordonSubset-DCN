clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
clim = [-0.5,0.5]; % color scale for FC

T = readtable('BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_stats_withMullen.csv');

load('BCP_Gordon_BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_7pt2min_randsample.mat')
zmatBCP = zmat;
for ii = 1:size(zmatBCP,1),zmatBCP(ii,ii,:) = 0;end
avg_zmatBCP = mean(zmatBCP,3);

load('washu120_parcellation_Gordon_20231101.mat')
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
%% try to do bootstrap CI for the average SI estimate for WU120
rng('default')
SilGordonMean = NaN(1000,1);
SilKardanMean = NaN(1000,1);
SilGordonSubsetMean =  NaN(1000,1);
SilKardanSubsetMean =  NaN(1000,1);

for isample = 1:1000
    isample
    idx = randsample(120,120,1);
    
    noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
    keepnets = IM.key(:,2)~=noneidx;
    M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    zmat_gordon_WashU120_subsample = mean(zmat_Gordon_all120(:,:,idx),3);
    D = calc_correlationdist(zmat_gordon_WashU120_subsample);
    s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
    SilGordonMean(isample) = mean(s);
    
    keepidx = SI_Infant_Gordon>0;
    s = silhouette_coef_mod(IM.key(keepidx,2),D(keepidx,keepidx),M);
    SilGordonSubsetMean(isample) = mean(s);
    
    noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
    keepnets = KardanIM.key(:,2)~=noneidx;
    M = ones(max(KardanIM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    zmat_Kardan_WashU120_subsample = mean(zmat_Kardan_all120(:,:,idx),3);
    D = calc_correlationdist(zmat_Kardan_WashU120_subsample);
    s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
    SilKardanMean(isample) = mean(s);
    
    keepidx = SI_Infant_Gordon>0;
    [~,sortid] = sort(IM.order);
    keepidx = keepidx(sortid);
    keepidx = keepidx(KardanIM.order);
    s = silhouette_coef_mod(KardanIM.key(keepidx,2),D(keepidx,keepidx),M);
    SilKardanSubsetMean(isample) = mean(s);
end
% save('./Results/WU120_bootstrap_results.mat','SilGordonMean','SilGordonSubsetMean','SilKardanMean','SilKardanSubsetMean')
%% Load WU120 bootstrap
% load('./Results/WU120_bootstrap_results.mat')

% figure;
% histogram(SilGordonMean-SilKardanMean);
% quantile(SilGordonMean-SilKardanMean,[0.025,0.975])
% quantile(SilGordonMean,[0.025,0.975])
% quantile(SilKardanMean,[0.025,0.975])
quantile(SilGordonSubsetMean,[0.025,0.975])
sum(SilGordonMean<0)
sum(SilKardanMean<0)
sum(SilGordonSubsetMean<0)
%% try to do bootstrap CI for the average SI estimate for BCP
rng('default')
SilGordonMean = NaN(1000,1);
SilKardanMean = NaN(1000,1);
SilGordonSubsetMean =  NaN(1000,1);
SilKardanSubsetMean =  NaN(1000,1);
SilGordon_BCP_eachROI = NaN(1000,333);
N = size(zmat_gordon_BCP,3)
for isample = 1:1000
    isample
    idx = randsample(N,N,1);
    
    zmat_gordon_BCP_subsample = mean(zmat_gordon_BCP(:,:,idx),3);
    D = calc_correlationdist(zmat_gordon_BCP_subsample);
    noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
    M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    
    keepidx = SI_Infant_Gordon>0;
    s = silhouette_coef_mod(IM.key(keepidx,2),D(keepidx,keepidx),M);
    SilGordonSubsetMean(isample) = mean(s);
    
    keepnets = IM.key(:,2)~=noneidx;
    s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
    SilGordonMean(isample) = mean(s);
    
    SilGordon_BCP_eachROI(isample,keepnets) = s;
    
    zmat_Kardan_BCP_subsample = mean(zmat_Kardan_BCP(:,:,idx),3);
    D = calc_correlationdist(zmat_Kardan_BCP_subsample);
    noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
    M = ones(max(KardanIM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    
    keepnets = KardanIM.key(:,2)~=noneidx;
    s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
    SilKardanMean(isample) = mean(s);
    
    keepidx = SI_Infant_Gordon>0;
    [~,sortid] = sort(IM.order);
    keepidx = keepidx(sortid);
    keepidx = keepidx(KardanIM.order);   
    s = silhouette_coef_mod(KardanIM.key(keepidx,2),D(keepidx,keepidx),M);
    SilKardanSubsetMean(isample) = mean(s);
end
% save('./results/BCP_bootstrap_results.mat','SilGordonMean','SilGordonSubsetMean','SilKardanMean','SilKardanSubsetMean','SilGordon_BCP_eachROI')
%% Load BCP bootstrap
% load('./results/BCP_bootstrap_results.mat')
% figure;
% histogram(SilGordonMean-SilKardanMean);
% quantile(SilGordonMean-SilKardanMean,[0.025,0.975])
% quantile(SilGordonMean,[0.025,0.975])
% quantile(SilKardanMean,[0.025,0.975])
quantile(SilKardanSubsetMean,[0.025,0.975])

sum(SilGordonMean<0)
sum(SilKardanMean<0)
sum(SilGordonSubsetMean<0)
%% Frequency of each ROI belonging the area subset with different bootstrapped sessions
% Supplementary Figure
freq_SIpos = mean(SilGordon_BCP_eachROI>0); % for each ROI, how many times in 1000 did it have SI>0. This is to see how robust our selected 166 area is

% plot frequency on the brain
[~,sortid] = sort(IM.order);
vals = freq_SIpos(sortid);
cmap = flipud(gray);

figure;
tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot_parcels_by_values(vals,'med',Parcels,[0,1],cmap) 
nexttile;
plot_parcels_by_values(vals,'lat',Parcels,[0,1],cmap) 
print(['./Figures/BrainPlotSIposfreq'],'-dtiff','-r300');
% Plot colorbar
[hCB,hf] = Plot.makecolorbar(cmap,[0,1],'h','frequency')
print(['./Figures/Colorbar'],'-dtiff','-r300');