%exported as 14" wide, 7" tall, 600dpi, fixed fonts at 7
clear all;

load mfb_workspace_titermatrix.mat;

%dctSitesNums = {'HC':1, 'KH':2, 'HU':3, 'DL':4, 'KG':5, 'ST':6, 'AG':7, 'DT':8, 'QN':9, 'BD':10 }
locnames = {'HC', 'KH', 'HU', 'DL', 'KG', 'ST', 'AG', 'DT', 'QN', 'BD' };

% --- COLUMNS of A ---
% 1     EXTRA01 - this is the local sample ID (i.e. last four digits) .. no need
%       to use this
% 2     EXTRA02 - this stores the location code; it's redundant
% 3     COLL_NUM
% 4     SAMPLE_ID,
% 5     AGE_MIN,
% 6     AGE_MAX,
% 7     GENDER_1ISMALE_VALUE,
% 8     DAY_VALUE,
% 9     MONTH_VALUE,
% 10    YEAR,
% 11    WARD_ID, 
% 12    EXTRA03, 
% 13    EXTRA04, 
% 14    EXTRA05,
% 15    SlideID,
% 16    BatchID,
% 17    SeriesID,
% 18    CollectionSite,
% 19    StudyCode,
% 20    SlideArrivalDay,
% 21    SlideArrivalMonth,
% 22    SlideArrivalYear,
% 23    ProcessedDay,
% 24    ProcessedMonth,
% 25    ProcessYear,
% 26    ScannedDay,
% 27    ScannedMonth,
% 28    ScannedYear,
% 29    SlideCorrectionFactor,
% 30    SlideAccepted,
% 31-34 :   H1_18_4F_TITER, H1_18_2F_TITER, H1_18_4F_TITERVALID, H1_18_2F_TITERVALID,
% 35-38 :   H1_77_4F_TITER, H1_77_2F_TITER, H1_77_4F_TITERVALID, H1_77_2F_TITERVALID,
% 39-42 :   H1_99_4F_TITER, H1_99_2F_TITER, H1_99_4F_TITERVALID, H1_99_2F_TITERVALID,
% 43-46 :   H1_07_4F_TITER, H1_07_2F_TITER, H1_07_4F_TITERVALID, H1_07_2F_TITERVALID,
% 47-50 :   H1_09_4F_TITER, H1_09_2F_TITER, H1_09_4F_TITERVALID, H1_09_2F_TITERVALID,
% 51-54 :   H3_68_4F_TITER, H3_68_2F_TITER, H3_68_4F_TITERVALID, H3_68_2F_TITERVALID,
% 55-58 :   H3_03_4F_TITER, H3_03_2F_TITER, H3_03_4F_TITERVALID, H3_03_2F_TITERVALID,
% 59-62 :   H3_05_4F_TITER, H3_05_2F_TITER, H3_05_4F_TITERVALID, H3_05_2F_TITERVALID,
% 63-66 :   H3_07_4F_TITER, H3_07_2F_TITER, H3_07_4F_TITERVALID, H3_07_2F_TITERVALID,
% 67-70 :   H3_09_4F_TITER, H3_09_2F_TITER, H3_09_4F_TITERVALID, H3_09_2F_TITERVALID,
% 71-74 :   H3_11_4F_TITER, H3_11_2F_TITER, H3_11_4F_TITERVALID, H3_11_2F_TITERVALID,
% 75-78 :   H3_13_4F_TITER, H3_13_2F_TITER, H3_13_4F_TITERVALID, H3_13_2F_TITERVALID,
% 79-82 :   H5_04_4F_TITER, H5_04_2F_TITER, H5_04_4F_TITERVALID, H5_04_2F_TITERVALID,
% 83-86 :   H5_07_4F_TITER, H5_07_2F_TITER, H5_07_4F_TITERVALID, H5_07_2F_TITERVALID,
% 87-90 :   H5_10_4F_TITER, H5_10_2F_TITER, H5_10_4F_TITERVALID, H5_10_2F_TITERVALID,
% 91-94 :   H7_03_4F_TITER, H7_03_2F_TITER, H7_03_4F_TITERVALID, H7_03_2F_TITERVALID,
% 95-98 :   H9_99_4F_TITER, H9_99_2F_TITER, H9_99_4F_TITERVALID, H9_99_2F_TITERVALID
% 99    :   AGE
% 100   :   FLOATING POINT SAMPLING DATE
% 101   :   FLOATING POINT BIRTHDAY


%%

for loc=1:10

    D  = SCORE( B(:,14)==loc, : );
    aa = B( B(:,14)==loc, 12 );
    aa = log2(aa);
    size(D,1);
    
    num_u5= size( aa( aa<log2(5) ), 1 );
    
    
    subplot(2,5,loc)
    plot( aa, D(:,1), 'o', 'MarkerSize',1 ); hold on;
    title(locnames(loc));
    text(log2(30),14,sprintf('n=%d',size(D,1)));
    text(log2(30),12,sprintf('u5=%d',num_u5));
    axis([-0.5 6.8 -15 15]);
    set(gca,'XTick', [log2(1) log2(2) log2(3) log2(4) log2(5) log2(10) log2(20) log2(40) log2(80) ] );
    set(gca,'XTickLabel', {'1','2','3','4','5','10','20','40','80'} );
    
    if loc==1 || loc==6
       ylabel('PC1 (11-dim PCA)'); 
    end

    
    [qq,ii]=sort(aa);
    aa_sorted=aa(ii);
    sc = D(:,1);
    sc_sorted = sc(ii);
    yy1 = smooth(aa_sorted,sc_sorted,0.5,'loess');
    plot(aa_sorted,yy1,'r-','LineWidth',2);

    axis([-2 log2(110) -1 25]);
    
end

%%

[COEFFH1,SCOREH1,latentH1,tsquaredH1,explainedH1,muH1] = pca(C(:,1:5), 'Centered', false);

[COEFFH3,SCOREH3,latentH3,tsquaredH3,explainedH3,muH3] = pca(C(:,6:11), 'Centered', false);


%%

AGE = B(:,12);
BIRTHDAY = B(:,13);
LOC = B(:,14);

save var_for_attackrate_fitting.mat AGE BIRTHDAY LOC COEFFH1 SCOREH1 latentH1 tsquaredH1 explainedH1 muH1 COEFFH3 SCOREH3 latentH3 tsquaredH3 explainedH3 muH3









