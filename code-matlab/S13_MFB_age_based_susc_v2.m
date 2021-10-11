clc;
clear all;
%==================

% 15*20 : fixed 36 half Dell monitor - 20pt H1N1, H3N2
% Illustrator: change H1N1, H3N2 to 26pt - Add locations

% MFB export:  7" by 9", fixed font at 10, custom renderer to vector

load var_for_attackrate_fitting.mat
locnames = {'HC', 'KH', 'HU', 'DL', 'KG', 'ST', 'AG', 'DT', 'QN', 'BD' };


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p,  [0.03 0.05], [0.05 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end


for idx = 1:8
    switch idx
        case 1
            SCORE = SCOREH1;
            lambda = 0.17392;
            lowerbound = 0.15926;
            upperbound = 0.18984;
            prms = [0.449821	7.74312	3.33754	0.173924	27658.9];
            ax = [-1 log2(85) -1.5 1.5];
            D  = SCORE( LOC==1 | LOC==2 | LOC==3, : );
            Age = AGE( LOC==1 | LOC==2 | LOC==3 );
            ytick = [-10:4:8];
            
        case 2
            SCORE = SCOREH3;
            lambda = 0.29519;
            lowerbound = 0.27609;
            upperbound = 0.31579;
            prms = [0.486546	11.5762	3.72952	0.295189	29460.1];
             ax = [-1 log2(85) -1 1.5];
            D  = SCORE( LOC==1 | LOC==2 | LOC==3, : );
            Age = AGE( LOC==1 | LOC==2 | LOC==3 );
            ytick = [-10:4:6];
        case 3
            SCORE = SCOREH1;
            lc = 1; % HCM
            lambda = 0.15352;
            lowerbound = 0.1307;
            upperbound = 0.16861;
            prms = [0.50321	7.6868	3.2778	0.15352	11291];
            ax = [-1 log2(85) -1.5 1.5];
            D  = SCORE( LOC==lc, : );
            Age = AGE( LOC==lc );
            ytick = [-6:4:8];
            
        case 4
            SCORE = SCOREH3;
            lc = 1; % HCM
            lambda = 0.29024;
            lowerbound = 0.26387;
            upperbound = 0.31874;
            prms = [0.721172	11.7831	3.90009	0.290242	12453.8];
             ax = [-1 log2(85) -1 1.5];
            D  = SCORE( LOC==lc, : );
            Age = AGE( LOC==lc );
            ytick = [-10:4:6];
            
        case 5
            SCORE = SCOREH1;
            lc = 2; %KH
            lambda = 0.14432;
            lowerbound = 0.11656;
            upperbound = 0.17872;
            prms = [0.47961	6.8694	3.3142	0.14432	8104];
            ax = [-1 log2(85) -1.5 1.5];
            D  = SCORE( LOC==lc, : );
            Age = AGE( LOC==lc );
            ytick = [-6:4:8];
            
        case 6
            SCORE = SCOREH3;
            lc = 2;%KH
            lambda = 0.28789;
            lowerbound = 0.25395;
            upperbound = 0.32615;
            prms = [0.202179	11.1418	3.51252	0.287893	8381.22];
             ax = [-1 log2(85) -1 1.5];
            D  = SCORE( LOC==lc, : );
            Age = AGE( LOC==lc );
            ytick = [-10:4:6];
            
        case 7
            SCORE = SCOREH1;
            lc = 3;%HU
            lambda = 0.24291;
            lowerbound = 0.20314;
            upperbound = 0.30197;
            prms = [0.40366	8.4113	3.4301	0.24791	8240.2];
            ax = [-1 log2(85) -1.5 1.5];
            D  = SCORE( LOC==lc, : );
            Age = AGE( LOC==lc );
            ytick = [-6:4:8];
            
        case 8
            SCORE = SCOREH3;
            lc = 3;%HU
            lambda = 0.32181;
            lowerbound = 0.27661;
            upperbound = 0.37469;
            prms = [0.480618	12.2934	3.67582	0.321811	8569.25];
             ax = [-1 log2(85) -1 1.5];
            D  = SCORE( LOC==lc, : );
            Age = AGE( LOC==lc );
            ytick = [-10:4:6];
    end
    
    ytick = [-10:4:8];
     
    PC1 = D(:,1);
    
    % this is a vector with ones and zeros
    idx_under_10 = Age<10;
    
    % the vector of ages under 10
    aa=Age(idx_under_10);

    % the vector of y observations under age 10
    %y_ob = log( PC1 - min(PC1) + 0.001 );
    y_ob = - log( max(PC1) - PC1 + 0.00001 );
    yy = y_ob(idx_under_10);
    
    [qq,ii]=sort(aa);
    aa_sorted=aa(ii);
    pc1_sorted = yy(ii);
    
    yy_sm = smooth(aa_sorted,pc1_sorted, 0.9 ,'loess');
    

    p=polyfit(aa_sorted,pc1_sorted,1);
    
    
    
    B = prms(1); % B is H
    C = prms(2); % C is K
    
    subplot(4,2,idx)
    %min(PC1);
    plot( aa, yy, 'LineStyle','none', 'Marker','o','MarkerEdgeColor',[.7 .7 .7], 'MarkerFaceColor', [.7 .7 .7],'MarkerSize',1); hold on;
    plot( aa_sorted, yy_sm, 'r-', 'LineWidth', 2);
    plot( aa_sorted, aa_sorted*p(1) + p(2), 'k-', 'LineWidth', 2);
    
    
    
%     x = [0.5:0.25:110];
%     y_ob = (B - PC1)/C
%     y_op = B - C*exp(-lambda.*x);
%     y_lo = B - C*exp(-lowerbound.*x);
%     y_up = B - C*exp(-upperbound.*x);
%     hold on;
%     log_x = log2(x);
%     
%     h = fill([log_x flip(log_x)], [y_lo flip(y_up)],[.4 .4 .4]);
%     alpha(0.5)
%     set(h,'edgecolor',[.4 .4 .4]);
%     
%     plot(log2(x),y_ob,'r','LineWidth',1.5)
      axis([0 10 -3 -1]);
    
%     set(gca, 'XTick', [log2(1) log2(2) log2(3) log2(4) log2(5) log2(10) log2(20) log2(40) log2(80)], 'XTickLabel', {'1','2','3','4','5','10','20','40','80'})
%     set(gca,'YTick',ytick)
    if idx == 7 || idx ==8
        xlabel('Age');
    end
   % ylabel('PC1')
    if idx == 1 || idx ==3|| idx ==5|| idx ==7
        ylabel(' log(max(PC1) - PC1) ');
    end
%     est = 1-exp(-lambda);
%     lb = 1-exp(-lowerbound);
%     ub = 1-exp(-upperbound);
%     
%     AR = ['AR: ' num2str(est*100,'%.1f') '%']
%     CI = ['95%CI: ' num2str(lb*100,'%.1f') '% - ' num2str(ub*100,'%.1f') '%']
%     text(log2(0.5)+0.3,6,AR)
%     text(log2(0.5)+0.3,4,CI)
    
    if idx == 1
        title('H1N1');
    end
    
    if idx == 2
        title('H3N2');
    end
    
   % title([num2str(est*100,'%.1f') '%: [' num2str(lb*100,'%.1f') '%, ' num2str(ub*100,'%.1f') '%]'])
end


