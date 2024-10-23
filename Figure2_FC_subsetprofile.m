clear;close all;clc;
addpath(genpath('./'));
%% ================ Figure 2 =======================
clim = [-0.5,0.5]; % color scale for FC
load(['/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/Parcels_','Gordon','.mat']);

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

%% Plot the subset of parcels on brain
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));

[~,sortid] = sort(IM.order);
vals = SI_Infant_Gordon(sortid);
rmidx = find(~(vals>=0));

load('MNI_coord_meshes_32k.mat');
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
Anat.CtxL.data=Parcel_Nets.CtxL; % plot Original Gordon Networks
Anat.CtxL.data(any(Parcels.CtxL==rmidx',2)) = 0;
Anat.CtxR.data=Parcel_Nets.CtxR;
Anat.CtxR.data(any(Parcels.CtxR==rmidx',2)) = 0;
params.Cmap.P=IM.cMap;
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'
params.lighting = 'gouraud';

figure;
tiledlayout(2,1,'TileSpacing','tight')
ax = nexttile;
params.view='lat';        % also, 'post','lat','med'
params.fig_handle = ax;
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);
ax = nexttile;
params.fig_handle = ax;
params.view='med';        % also, 'post','lat','med'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);
% print(['./Figures/GordonNetwork'],'-dtiff','-r300');
%% Get number of nodes and average FC for before and after removal
Nnet = 13;
Nsubj = size(zmatBCP,3)
[within_net_FC,within_net_FC_partial,within_net_FC_excluded] = deal(NaN(Nnet,Nsubj));
[count_node,count_node_partial] = deal(NaN(Nnet,1));
for inet = setdiff(1:Nnet,noneidx)
    inet
    netidx = IM.key(:,2)==inet;
    netidx_partial = IM.key(:,2)==inet & SI_Infant_Gordon>0;
    netidx_excluded = IM.key(:,2)==inet & SI_Infant_Gordon<0;
    count_node(inet) = sum(netidx);
    count_node_partial(inet) = sum(netidx_partial);
    for isubj = 1:Nsubj
        UDidx = get_triu_idx(sum(netidx));
        tmp = zmat_gordon_BCP(netidx,netidx,isubj);
        within_net_FC(inet,isubj) = mean(tmp(UDidx));
        
        UDidx = get_triu_idx(sum(netidx_partial));
        tmp = zmat_gordon_BCP(netidx_partial,netidx_partial,isubj);
        within_net_FC_partial(inet,isubj) = mean(tmp(UDidx));  
        
        UDidx = get_triu_idx(sum(netidx_excluded));
        tmp = zmat_gordon_BCP(netidx_excluded,netidx_excluded,isubj);
        within_net_FC_excluded(inet,isubj) = mean(tmp(UDidx));   
    end
end
%% number of nodes remained in each network
figure('units','inches','position',[1,1,6,3]);hold on;
icounter = 0;
clear hb
for inet = setdiff(1:Nnet,noneidx)
    icounter = icounter+1;
    hb(1) = bar(icounter-0.2,count_node(inet),0.2,'FaceColor',IM.cMap(inet,:),'EdgeColor',IM.cMap(inet,:));
    hb(2) = bar(icounter+0.2,count_node_partial(inet),0.2,'FaceColor','w','EdgeColor',IM.cMap(inet,:));
end
xticks(1:Nnet-1);xticklabels(IM.Nets(setdiff(1:Nnet,noneidx)));xtickangle(30);
set(gca,'FontWeight','Bold','FontSize',10);
grid minor
legend(hb,{'All','Subset'},'location','eastoutside');
ylabel('Number of parcels');
% print('./Figures/Gordon_BCP_numnodes_original_partial_barchart','-dpdf');

figure;
subplot(1,2,1);
cMap  = IM.cMap(count_node>0,:);
h = pie(count_node(count_node>0));
for i = 1:size(cMap,1)
    h(2*(i-1)+1).FaceColor = cMap(i,:);
end
title('All');
set(gca,'FontSize',12);
subplot(1,2,2);
cMap  = IM.cMap(count_node_partial>0,:);
h = pie(count_node_partial(count_node_partial>0));
for i = 1:size(cMap,1)
    h(2*(i-1)+1).FaceColor = cMap(i,:);
end
title('Subset');
set(gca,'FontSize',12);
% print('./Figures/Gordon_BCP_numnodes_original_partial_piechart','-dpdf');
%% Boxplot of within network FC in All VS Subset (Infant)
retained_nets = [1,2,4,5,6,9,12,13]
clear h
[~,p]=arrayfun(@(i)ttest(within_net_FC(i,:)',within_net_FC_partial(i,:)'),retained_nets);
p = mafdr(p,'BHFDR',true);
figure('units','inches','position',[1,1,7,3]);hold on;
X = [within_net_FC(retained_nets,:)',within_net_FC_partial(retained_nets,:)'];
for k = 1:length(retained_nets)
    plot([k-0.2,k+0.2],X(:,[k,k+8])','Color',[0.5,0.5,0.5,0.1])
end
boxplot(X,1:16,'position',[[1:8]-0.2,[1:8]+0.2],'colors',repelem([IM.cMap(retained_nets,:)],2,1),'symbol','.'); 
h = findobj(gca, 'Tag', 'Box');
h2 = findobj(gca,'Tag','Upper Whisker');
h3 = findobj(gca,'Tag','Lower Whisker');
for k =  2:2:length(h)
    h(k).LineWidth = 2;
    h2(k).LineWidth = 2;
    h3(k).LineWidth = 2;
end

sigstar(arrayfun(@(i)[i-0.2,i+0.2],1:length(retained_nets),'UniformOutput',false),p)

xticks(1:length(retained_nets))
xticklabels(IM.Nets(retained_nets))
xtickangle(30);
set(gca,'FontSize',12);
title('Within network FC');
% grid minor
ylabel('z(r)')
% print('./Figures/Gordon_BCP_withinnetFC_original_partial','-dtiff','-r300');

%% Is this difference between within network using All V.S. Subset related to age?
[r,p] = deal(NaN(size(retained_nets)));
counter = 0;
for inet = retained_nets
    counter = counter+1;
    figure('units','inches','position',[1,1,3,3]);hold on;
    scatter(T.age_yrs,within_net_FC(inet,:)-within_net_FC_partial(inet,:),15,IM.cMap(inet,:),'.');
    [r(counter),p(counter)] = corr(T.age_yrs,[within_net_FC(inet,:)-within_net_FC_partial(inet,:)]');
    h = refline; h.Color = 'r';
    xlabel('Age (yrs)');
    set(gca,'FontSize',12);
    text(0.1,0.1,sprintf('r = %1.2f, p = %1.3f',r(counter),p(counter)),'Units','Normalized')
    title(IM.Nets{inet});
    print(['./Figures/Individual_withinFC_Kardan_Infant_',IM.Nets{inet},'.tif'],'-dtiff','-r300');
end
close all
FDR = mafdr(p,'BHFDR',true)