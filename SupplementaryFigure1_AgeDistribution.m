%% Get BCP age
addpath(genpath('./'))
T = readtable('BCP_Jan2023_QCpass_asleep_atleast7pt2min_UNC_UMN_20240124_stats_withMullen.csv');
%%
[sub_id,b] = unique(T.sub_id,'stable');
[~,sortid] = sort(T.sex_M1F0(b));
figure('position',[100 100 600 400]);hold on;
counter = 0;
agediff = [];
longitudinal_points = [];
for i = sortid'
    counter = counter+1;
    age_yrs = sort(T.age_yrs(string(T.sub_id)==sub_id{i}));
    longitudinal_points(i) = length(age_yrs);
    if length(age_yrs)>1
        agediff = [agediff;diff(age_yrs)];
    end
    if unique(T.sex_M1F0(string(T.sub_id)==sub_id{i}))
        h(1) = plot(age_yrs,repelem(counter,length(age_yrs)),'.-','Color','b');
    else
        h(2)= plot(age_yrs,repelem(counter,length(age_yrs)),'.-','Color','r');
    end
end
legend(h,{'Male','Female'},'location','eastoutside');
xlabel('yrs');
ylabel('subjects');
title('BCP');
set(gca,'FontSize',15);
ylim([0,length(sub_id)+1]);
%% Longitudinal points
[num,lab] = histcounts(categorical(longitudinal_points));
figure('position',[100 100 600 400]);hold on;
bar(1:length(num),num);
for i = 1:length(num)
    text(i-0.1,num(i)+3,string(num(i)),'FontSize',15);
end
xticks(1:length(num))
xticklabels(lab);
ylim([0,95])
set(gca,'FontSize',15);
xlabel('Number of longitudinal points')
ylabel('Count');