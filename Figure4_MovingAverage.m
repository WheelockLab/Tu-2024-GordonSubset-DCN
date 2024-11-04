clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
clim =[-0.5,0.5];
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
%% Moving average window for age

noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
keepnets = IM.key(:,2)~=noneidx;
M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));

[Ages,AgesortID] = sort(T.age_yrs);
Ages = Ages(1:281);AgesortID = AgesortID(1:281); % after 281 the data point is very scattered
windowsz = 20;% 20 sessions
windowstep = 1;
N = floor((length(Ages)-windowsz)/windowstep);
[sil_Gordon,sil_Gordon_subset,sil_Kardan,sil_Kardan_subset,agemean] = deal(NaN(N,1));
[sil_Gordon_bstrp,sil_Gordon_subset_bstrp,sil_Kardan_bstrp,sil_Kardan_subset_bstrp] = deal(NaN(N,1000));

for jj =1:N
    jj
    agemean(jj) =mean(Ages([1:windowsz]+(jj-1)*windowstep)); % ages weighted by the mean
    
    wholesample = AgesortID([1:windowsz]+(jj-1)*windowstep);
    tmpzmat=mean(zmat_gordon_BCP(:,:,wholesample),3);
    D = calc_correlationdist(tmpzmat);
    s = silhouette_coef_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
    sil_Gordon(jj) = mean(s);
    SI_Adult_Gordon_tmp = NaN(Nroi,1);
    SI_Adult_Gordon_tmp(keepnets) = s;
    s = silhouette_coef_mod(IM.key(SI_Infant_Gordon>0,2),D(SI_Infant_Gordon>0,SI_Infant_Gordon>0),M);
    sil_Gordon_subset(jj) = mean(s);% mean(SI_Adult_Gordon_tmp(SI_Infant_Gordon>0)); %20240509: recalculate with the remaining parcels instead
    
end

% repeat for Kardan ?
noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
keepnets = KardanIM.key(:,2)~=noneidx;
M = ones(max(KardanIM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));

[~,sortid] = sort(IM.order);
keepidx = SI_Infant_Gordon>0;
keepidx = keepidx(sortid);
keepidx = keepidx(KardanIM.order);

for jj =1:N
    jj
    wholesample = AgesortID([1:windowsz]+(jj-1)*windowstep);
    tmpzmat=mean(zmat_Kardan_BCP(:,:,wholesample),3);
    D = calc_correlationdist(tmpzmat);
    s = silhouette_coef_mod(KardanIM.key(keepnets,2),D(keepnets,keepnets),M);
    sil_Kardan(jj) = mean(s);
    s = silhouette_coef_mod(KardanIM.key(keepidx,2),D(keepidx,keepidx),M);
    sil_Kardan_subset(jj) = mean(s);% mean(SI_Adult_Gordon_tmp(SI_Infant_Gordon>0)); %20240509: recalculate with the remaining parcels instead
end

save('./Results/moving_avg_results.mat','sil_Gordon*','sil_Kardan*','agemean')
%% Plot moving average (Figure 4A)
legendstr = {'Gordon','Gordon (Subset)','Kardan','Kardan (Subset)'}
clear h
figure('position',[100 100 400 400]);hold on;
h(1) = plot(agemean,sil_Gordon,'LineWidth',2);
h(2) = plot(agemean,sil_Gordon_subset,'LineWidth',2);
h(3) = plot(agemean,sil_Kardan,'LineWidth',2);
h(4) = plot(agemean,sil_Kardan_subset,'LineWidth',2);
xlim([0.5,2.5]);
Plot.vline(agemean(53)); %~ 1 year
Plot.vline(agemean(223)); % ~2  year
legend(h,legendstr,'location','best');
legend('location','southoutside');
ylabel('SI');
xlabel('Age (yrs)')
set(gca,'FontSize',12);
print(gcf,'./Figures/SilMovingAvg','-dpdf'); 
% print(gcf,'./Figures/SilMovingAvg','-dtiff','-r300'); 

%% Plot examples for 1 year and 2 year (Figure 4B)
keepidx = SI_Infant_Gordon>0;
[~,sortid] = sort(IM.order);
keepidx2 = keepidx(sortid);
keepidx2 = keepidx2(KardanIM.order);
for jj = [53,223]
    figure('position',[100 100 800 200])
    noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
    keepnets = IM.key(:,2)~=noneidx;
    
    tmpzmat=mean(zmat_gordon_BCP(:,:,AgesortID([1:windowsz]+(jj-1)*windowstep)),3);

    subplot(1,4,1);
    Matrix_Org3(tmpzmat(keepnets,keepnets),...
    IM.key(keepnets,:),10,clim,IM.cMap,0);
    
    text(0.2,-0.1,sprintf('avg SI = %2.3f',sil_Gordon(jj)),'Units','normalized','FontWeight','Bold','FontSize',12);
    
    subplot(1,4,2);
    
    Matrix_Org3(tmpzmat(keepidx,keepidx),...
    IM.key(keepidx,:),10,clim,IM.cMap,0);
    text(0.2,-0.1,sprintf('avg SI = %2.3f',sil_Gordon_subset(jj)),'Units','normalized','FontWeight','Bold','FontSize',12);
    
    noneidx = find((string(KardanIM.Nets)=="None")|(string(KardanIM.Nets)=="USp"));
    keepnets = KardanIM.key(:,2)~=noneidx;
    tmpzmat=mean(zmat_Kardan_BCP(:,:,AgesortID([1:windowsz]+(jj-1)*windowstep)),3);
    subplot(1,4,3);
    Matrix_Org3(tmpzmat(keepnets,keepnets),...
        KardanIM.key(keepnets,:),10,clim,KardanIM.cMap,0);
    text(0.2,-0.1,sprintf('avg SI = %2.3f',sil_Kardan(jj)),'Units','normalized','FontWeight','Bold','FontSize',12);
    
    subplot(1,4,4);
    Matrix_Org3(tmpzmat(keepidx2,keepidx2),...
        KardanIM.key(keepidx2,:),10,clim,KardanIM.cMap,0);
    text(0.2,-0.1,sprintf('avg SI = %2.3f',sil_Kardan_subset(jj)),'Units','normalized','FontWeight','Bold','FontSize',12);
    print(['./Figures/DifferentSorting_',num2str(jj)],'-dpdf'); % uncomment relevant sections to also plot the spheres
end