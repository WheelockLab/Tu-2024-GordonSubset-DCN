clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
load('BCP_Gordon_BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_7pt2min_randsample.mat')

zmatBCP = zmat;
for ii = 1:size(zmatBCP,1),zmatBCP(ii,ii,:) = 0;end
avg_zmatBCP = mean(zmatBCP,3);

load('IM_Gordon_13nets_333Parcels_renamed.mat','IM')
zmat_gordon_BCP = zmatBCP(IM.order,IM.order,:);
avg_zmat_gordon_BCP = avg_zmatBCP(IM.order,IM.order);
[~,Nroi]=size(zmat_gordon_BCP);

% load Dworestky ROIs
Dworetsky153ROI = readtable('/data/wheelock/data1/software/Dworetsky_etal_ConsensusNetworks/153ProbabilisticROIs_MNI_info.txt');% Dworetsky et al. 153 high probability ROI across individuals
%% Get SI for each parcel based on Gordon parcellations on BCP data
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
D = calc_correlationdist(avg_zmat_gordon_BCP);
s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);

SI_Infant_Gordon = NaN(Nroi,1);
SI_Infant_Gordon(keepnets) = s;
%% To-do: add in the plot on brain surface part for Figure 5A/B

%% ========== Figure 5C ==========
%% Statistically compare how similar they are to Dworestky ROIs
[~,sortid] = sort(IM.order);
GordonROIcoord = IM.ROIxyz(sortid,:);
vals = SI_Infant_Gordon(sortid);
rmidx = find(isnan(vals)|vals<0);
rmidx2 = find(vals<0); % not including the None networks
keepidx = setdiff(1:333,rmidx);


D = pdist2(Dworetsky153ROI{:,2:4},GordonROIcoord(keepidx,:),'euclidean');
realD = min(D,[],2);
D = pdist2(Dworetsky153ROI{:,2:4},GordonROIcoord(rmidx2,:),'euclidean');
realDrm = min(D,[],2);
realDiff = mean(realDrm-realD);

rng('default');

nonnanidx = find(~isnan(vals));
[randDiff] = deal(NaN(length(realDiff),1000));
[randDrm_all,randD_all] = deal(NaN(length(realD),1000));
for iter = 1:1000
    rrmidx = randsample(nonnanidx,length(rmidx2),0);  
    rkeepidx = setdiff(nonnanidx,rrmidx);
    D = pdist2(Dworetsky153ROI{:,2:4},GordonROIcoord(rkeepidx,:),'euclidean');
    randD= min(D,[],2);
    D = pdist2(Dworetsky153ROI{:,2:4},GordonROIcoord(rrmidx,:),'euclidean');
    randDrm = min(D,[],2);
    randDiff(iter) = mean(randDrm-randD);
    randDrm_all(:,iter) = randDrm;
    randD_all(:,iter) = randD;
end

figure('position',[100 100 150,250]);
histogram(randDiff);
xline(realDiff,'LineStyle','--');
yticks([]);
% title('Average Euclidean Distance to Dworetsky ROIs');
xlabel('mm');
set(gca,'FontSize',12,'FontWeight','Bold');
% print('./Figures/EuclideanDistanceToDworetskyROI','-dpdf')
%% Real (Supplementary Figure)
figure;hold on;
boxplot([realD,realDrm]);
plot([realD,realDrm]','Color',[0.5 0.5 0.5 0.3])
xticklabels({'Subset','Alternative'});
ylabel('Distance (mm)')
set(gca,'FontSize',12);
print('./Figures/EuclideanDistanceToDworetskyROI_real','-dpng')
%% Example Null
% figure;hold on;
% histogram(randD_all(:,iter),0:5:35,'EdgeColor','None');
% histogram(randDrm_all(:,iter),0:5:35,'EdgeColor','None');
iter = 2
figure;hold on;
boxplot([randD_all(:,iter),randDrm_all(:,iter)]);
plot([randD_all(:,iter),randDrm_all(:,iter)]','Color',[0.5 0.5 0.5 0.3])
xticklabels({'Subset','Alternative'});
ylabel('Distance (mm)')
set(gca,'FontSize',12);

%% Plot individual distributions
figure;
subplot(1,2,1);
histogram(randD);
xline(mean(realD),'LineStyle','--');
subplot(1,2,2);
histogram(randDrm);
xline(mean(realDrm),'LineStyle','--');

%% Save coordinates to csv
writematrix(GordonROIcoord(keepidx,:),'./GordonSubsetCoordinates.csv')
writematrix(keepidx','./GordonSubsetIdxIn333.txt')