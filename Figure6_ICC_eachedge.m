clear;close all;clc;
addpath(genpath('./'));
%% Set some parameters and load data
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

%% Reliability for each network
load('./Results/zmat_Gordon_unsorted_BCP_6min.mat','zmat1','zmat2','rmidx')
Npar = 333;
UDidx = get_triu_idx(Npar);
[~,sortid] = sort(IM.order);
key = IM.key(:,2); 

zmat1_square =zeros(Npar,Npar,size(zmat1,2));
for ii = 1:size(zmat1,2)
    tmp = zeros(Npar,Npar);
    tmp(UDidx) = zmat1(:,ii);
    tmp = tmp+tmp';
    zmat1_square(:,:,ii) = tmp;
end
zmat1_square_sorted = zmat1_square(IM.order,IM.order,:);

zmat2_square =zeros(Npar,Npar,size(zmat2,2));
for ii = 1:size(zmat2,2)
    tmp = zeros(Npar,Npar);
    tmp(UDidx) = zmat2(:,ii);
    tmp = tmp+tmp';
    zmat2_square(:,:,ii) = tmp;
end
zmat2_square_sorted = zmat2_square(IM.order,IM.order,:);

[allicc,alllb,allub,allp] = deal(NaN(max(IM.key(:,2)),1));
for inet = setdiff(unique(key),0)'
    data = reshape(zmat1_square_sorted(key==inet,key==inet,:),[],size(zmat1_square_sorted,3));
    N =sum(data~=0);
    data = [sum(data)./N]';
    data2 = reshape(zmat2_square_sorted(key==inet,key==inet,:),[],size(zmat2_square_sorted,3));
    data2 = [sum(data2)./N]';
    
     [r, LB, UB, F, ~,~, p] = ICC([data,data2], 'C-1');
     allicc(inet)=r;
     alllb(inet)=LB;
     allub(inet)=UB;
     allp(inet)=p;
end

[~,sortid] = sort(IM.order);
keysubset = IM.key(:,2); keysubset(SI_Infant_Gordon<0|isnan(SI_Infant_Gordon)) = 0;

[allicc2,alllb2,allub2,allp2] = deal(NaN(max(keysubset),1));
for inet = setdiff(unique(keysubset),0)'
        data = reshape(zmat1_square_sorted(keysubset==inet,keysubset==inet,:),[],size(zmat1_square_sorted,3));
        N =sum(data~=0);
        data = [sum(data)./N]';
        data2 = reshape(zmat2_square_sorted(keysubset==inet,keysubset==inet,:),[],size(zmat2_square_sorted,3));
        data2 = [sum(data2)./N]';
        [r, LB, UB, F, ~,~, p] = ICC([data,data2], 'C-1');
        allicc2(inet)=r;
        alllb2(inet)=LB;
        allub2(inet)=UB;
        allp2(inet)=p;   
end
%% ICC of each edge
r = zeros(size(zmat1_square_sorted,1));
for ii = 1:size(zmat1_square_sorted,1)
    for jj = ii+1:size(zmat1_square_sorted,1)
        data = reshape(zmat1_square_sorted(ii,jj,:),size(zmat1_square_sorted,3),[]);
        data2 = reshape(zmat2_square_sorted(ii,jj,:),size(zmat2_square_sorted,3),[]);
        r(ii,jj)= ICC([data,data2], 'C-1');
    end
end
r = r+r';

netidx = [1,2,4,5,6,9,12,13];%setdiff(unique(keysubset),0)';
counter = 0;
[m1,m2,s1,s2,p,d] = deal(NaN(max(key),1));
for inet = netidx
    counter = counter+1;
    x1 = nonzeros(r(key==inet,key==inet));
    x2 = nonzeros(r(keysubset==inet,keysubset==inet));
    m1(inet) = mean(x1);
    m2(inet) = mean(x2);
    s1(inet) = std(x1)/sqrt(length(x1));
    s2(inet) = std(x2)/sqrt(length(x2));
    [~,p(inet)] = ttest2(x1,x2);
    d(inet) = computeCohen_d(x1, x2);
end

p = mafdr(p,'BHFDR',true);

figure('units','inches','position',[1,1,7,3]);hold on;
icounter = 0;
clear hb
cumpairs = {};cump = [];
for inet = netidx%setdiff(1:Nnet,noneidx)
    icounter = icounter+1;
    hb(1) = bar(icounter-0.2,m1(inet),0.2,'FaceColor',IM.cMap(inet,:),'EdgeColor',IM.cMap(inet,:));
    hb(2) = bar(icounter+0.2,m2(inet),0.2,'FaceColor','w','EdgeColor',IM.cMap(inet,:));
    errorbar(icounter-0.2,m1(inet),s1(inet),'Color','k','CapSize',0);
    errorbar(icounter+0.2,m2(inet),s2(inet),'Color','k','CapSize',0);
    cumpairs{end+1} = [icounter-0.2,icounter+0.2];
    cump(end+1) = p(inet);
end
sigstar(cumpairs,cump)
xticks(1:length(netidx));xticklabels(IM.Nets(netidx));xtickangle(30);
set(gca,'FontWeight','Bold','FontSize',10);
grid minor
legend(hb,{'All','Subset'},'location','eastoutside');
ylabel('ICC of within network edges');
print(gcf,'./Figures/withinnetworkedgeICC','-dtiff','-r300');