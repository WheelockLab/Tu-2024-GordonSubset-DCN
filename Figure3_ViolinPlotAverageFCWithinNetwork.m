% Figure3_ViolinPlotAverageFCWithinNetwork
clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
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

[~,Nroi,Nsubj]=size(zmat_gordon_BCP);
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
%% Get number of nodes and average FC for before and after removal
Nnet = 13;
[within_net_FC,within_net_FC_partial,within_net_FC_excluded] = deal(NaN(Nnet,Nsubj));
[within_net_FC_Adult,within_net_FC_Adult_partial] = deal(NaN(Nnet,size(zmat_gordon_WashU120,3)));
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
for inet = setdiff(1:Nnet,noneidx)
    inet
    netidx = IM.key(:,2)==inet;
    netidx_partial = IM.key(:,2)==inet & SI_Infant_Gordon>0;
    for isubj = 1:size(zmat_Gordon_all120,3)
        UDidx = get_triu_idx(sum(netidx));
        tmp = zmat_Gordon_all120(netidx,netidx,isubj);
        within_net_FC_Adult(inet,isubj) = mean(tmp(UDidx));
        
        UDidx = get_triu_idx(sum(netidx_partial));
        tmp = zmat_Gordon_all120(netidx_partial,netidx_partial,isubj);
        within_net_FC_Adult_partial(inet,isubj) = mean(tmp(UDidx));
    end
end
%% Compare withinnet adult and infant
clear h
currnets = [1,2,4,5,6,9,12,13]
[~,p]=arrayfun(@(i)ttest2(within_net_FC_Adult(i,:)',within_net_FC(i,:)'),currnets);
p = mafdr(p,'BHFDR',true);
figure('units','inches','position',[1,1,5.5,3]);hold on;
icounter = 0;

for inet = currnets
    icounter = icounter+1;
    hh =Violin(within_net_FC_Adult(inet,:), icounter -0.2,'ViolinColor',IM.cMap(inet,:),'Width',0.1);
    h(1) = hh.ViolinPlot;
end
icounter = 0;
for inet =currnets
    icounter = icounter+1;
    try
        hh = Violin(within_net_FC(inet,:), icounter +0.2,'ViolinColor',[0.8 0.8 0.8],'EdgeColor',IM.cMap(inet,:),'Width',0.1);
        h(2)  = hh.ViolinPlot;
    catch
    end
end
xticks(1:Nnet-1);xticklabels(IM.Nets(currnets ));xtickangle(30);
set(gca,'FontWeight','Bold','FontSize',10);
title('Within-network FC');
grid minor
xlim([0,length(currnets)+1])
ylabel('z(r)')
sigstar(arrayfun(@(i)[i-0.2,i+0.2],1:length(p),'UniformOutput',false),p)
% legend(h,{'Adult (All)','Infant (All)'},'location','eastoutside')
ylim([0,1])

print('./Figures/Gordon_BCP_withinnetFC_adult_infant_all','-dpdf');

d = NaN(13,1);
for inet = currnets
    d(inet) = computeCohen_d(within_net_FC_Adult(inet,:),within_net_FC(inet,:),'independent');
end
%% Compare withinnet adult and infant (subset)
clear h
currnets = [1,2,4,5,6,9,12,13]
[~,p]=arrayfun(@(i)ttest2(within_net_FC_Adult_partial(i,:)',within_net_FC_partial(i,:)'),currnets);

p = mafdr(p,'BHFDR',true);
figure('units','inches','position',[1,1,5.5,3]);hold on;
icounter = 0;

for inet = currnets
    icounter = icounter+1;
    hh =Violin(within_net_FC_Adult_partial(inet,:), icounter -0.2,'ViolinColor',IM.cMap(inet,:),'Width',0.1);
    h(1) = hh.ViolinPlot;
end
icounter = 0;
for inet =currnets
    icounter = icounter+1;
    try
        hh = Violin(within_net_FC_partial(inet,:), icounter +0.2,'ViolinColor',[0.8 0.8 0.8],'EdgeColor',IM.cMap(inet,:),'Width',0.1);
        h(2)  = hh.ViolinPlot;
    catch
    end
end
xticks(1:Nnet-1);xticklabels(IM.Nets(currnets ));xtickangle(30);
set(gca,'FontWeight','Bold','FontSize',10);
title('Within-network FC');
grid minor
xlim([0,length(currnets)+1])
ylabel('z(r)')
sigstar(arrayfun(@(i)[i-0.2,i+0.2],1:length(p),'UniformOutput',false),p)
% legend(h,{'Adult (Subset)','Infant (Subset)'},'location','eastoutside')
ylim([0,1])

print('./Figures/Gordon_BCP_withinnetFC_adult_infant_subset','-dpdf');

d = NaN(13,1);
for inet = currnets
    d(inet) = computeCohen_d(within_net_FC_Adult_partial(inet,:),within_net_FC_partial(inet,:),'independent');
end