clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
load('BCP_Gordon_BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_7pt2min_randsample.mat')

zmatBCP = zmat;
for ii = 1:size(zmatBCP,1),zmatBCP(ii,ii,:) = 0;end
avg_zmatBCP = mean(zmatBCP,3);

load('IM_Gordon_13nets_333Parcels_renamed.mat','IM','Parcel_Nets')
zmat_gordon_BCP = zmatBCP(IM.order,IM.order,:);
avg_zmat_gordon_BCP = avg_zmatBCP(IM.order,IM.order);
[~,Nroi]=size(zmat_gordon_BCP);

% load Dworestky ROIs
Dworetsky153ROI = readtable('/data/wheelock/data1/software/Dworetsky_etal_ConsensusNetworks/153ProbabilisticROIs_MNI_info.txt');% Dworetsky et al. 153 high probability ROI across individuals

% Load parcel positions
load(['Parcels_','Gordon','.mat']);
%% Get SI for each parcel based on Gordon parcellations on BCP data
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
D = calc_correlationdist(avg_zmat_gordon_BCP);
s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);

SI_Infant_Gordon = NaN(Nroi,1);
SI_Infant_Gordon(keepnets) = s;
%% Plot the subset of parcels on brain
[~,sortid] = sort(IM.order);
vals = SI_Infant_Gordon(sortid);

% removeidx = find(isnan(vals)|vals<0); % Only Keep Area Subset
removeidx = find(isnan(vals)|vals>0);% Only Keep Alternative Areas

% rmidx = NaN;
[Lindtrunc] = with_without_mw_conversion('Lindtrunc');
[Rindtrunc] = with_without_mw_conversion('Rindtrunc');

load('MNI_coord_meshes_32k.mat');
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
Anat.CtxL.data=Parcel_Nets.CtxL; % plot Original Gordon Networks
Anat.CtxL.data(any(Parcels.CtxL==removeidx',2)) = 0; % get rid of the keepidx
Anat.CtxR.data=Parcel_Nets.CtxR;
Anat.CtxR.data(any(Parcels.CtxR==removeidx',2)) = 0; % get rid of the keepidx
params.Cmap.P=IM.cMap;
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'
params.lighting = 'gouraud';
 
Anat2.CtxL=MNIl;Anat2.CtxR=MNIr;
Anat2.ctx = 'inf';
ROI.coord = Dworetsky153ROI{:,2:4};
ROI.radius = repmat(3.5,153,1); 
ROI.color = repmat([188 143 143]/255,153,1);
ROI.Network = Dworetsky153ROI{:,6};

figure;
tiledlayout(2,1,'tilespacing','tight')
ax = nexttile;%subplot(2,1,1);
params.view='lat';        % also, 'post','lat','med'
params.fig_handle = ax;
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);

% draw spheres
Anat2.view ='lat';
foci.lighting=params.lighting;
for j=1:length(unique(ROI.Network))
    keep=find(ROI.Network==j);
    foci.color=ROI.color(keep,:);
    foci.radius=ROI.radius(keep,:);
    foci.location=ROI.coord(keep,:);
    foci=AdjustFoci(foci,Anat2.CtxL,Anat2.CtxR,Anat2.view,Anat2.ctx);
    hold on
    Draw_Foci(foci,10)
end

% ax = subplot(2,1,2);
ax = nexttile;
params.fig_handle = ax;
params.view='med';        % also, 'post','lat','med'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);

% draw spheres
Anat2.view ='med';
foci.lighting=params.lighting;
for j=1:length(unique(ROI.Network))
    keep=find(ROI.Network==j);
    foci.color=ROI.color(keep,:);
    foci.radius=ROI.radius(keep,:);
    foci.location=ROI.coord(keep,:);
    foci=AdjustFoci(foci,Anat2.CtxL,Anat2.CtxR,Anat2.view,Anat2.ctx);
    hold on
    Draw_Foci(foci,10)
end
% draw spheres
Anat2.view ='med';
foci.lighting=params.lighting;
for j=1:length(unique(ROI.Network))
    keep=find(ROI.Network==j);
    foci.color=ROI.color(keep,:);
    foci.radius=ROI.radius(keep,:);
    foci.location=ROI.coord(keep,:);
    foci=AdjustFoci(foci,Anat2.CtxL,Anat2.CtxR,Anat2.view,Anat2.ctx);
    hold on
    Draw_Foci(foci,10)
end
print('./Figures/GordonSubset','-dtiff','-r300')

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
figure('position',[100 100 300 400]);hold on;
boxplot([realD,realDrm],'plotstyle','compact','labels',{'Subset','Alternative'});
% plot([realD,realDrm]','Color',[0.5 0.5 0.5 0.3])
% xticklabels({'Subset','Alternative'});
ylabel('Distance (mm)')
set(gca,'FontSize',12);
print('./Figures/EuclideanDistanceToDworetskyROI_real','-dpng','-r300')
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