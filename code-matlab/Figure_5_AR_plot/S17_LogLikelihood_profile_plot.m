clc;
clear all;
%==============================================

load ../var_for_attackrate_fitting.mat
locnames = {'HC', 'KH', 'HU', 'DL', 'KG', 'ST', 'AG', 'DT', 'QN', 'BD' };


subtype = 'H1N1';
%subtype = 'H3N2';


panel = 1;
for lc=1:3
    llh_profile = dlmread(['./' subtype '/LLHProfile_' locnames{lc} '.txt'])
    subplot(2,3,panel)
    plot(llh_profile(:,4),llh_profile(:,5),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r');
    hold on;
    plot([min(llh_profile(:,4))-0.01 max(llh_profile(:,4))+0.01],[min(llh_profile(:,5))+1.92 min(llh_profile(:,5))+1.92])
    
    [est, lb, ub] = Find95CI(llh_profile(:,4),llh_profile(:,5))
    
    title([locnames{lc} ' - ' num2str(est) ': [' num2str(lb) ', ' num2str(ub) ']'])
    panel = panel + 1;
end



subtype = 'H3N2';
for lc=1:3
    llh_profile = dlmread(['./' subtype '/LLHProfile_' locnames{lc} '.txt'])
    subplot(2,3,panel)
    plot(llh_profile(:,4),llh_profile(:,5),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r');
    hold on;
    plot([min(llh_profile(:,4))-0.01 max(llh_profile(:,4))+0.01],[min(llh_profile(:,5))+1.92 min(llh_profile(:,5))+1.92])
    
    [est, lb, ub] = Find95CI(llh_profile(:,4),llh_profile(:,5))
    
    title([locnames{lc} ' - ' num2str(est) ': [' num2str(lb) ', ' num2str(ub) ']'])
    panel = panel + 1;
end