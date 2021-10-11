clear all;

global PC1 Age;

B = csvread('../data/TITERS_02FL_CompleteSetOf11HumanAg_Compact.20201120.csv',1,0);
%B = csvread('../data/TITERS_02FL_CompleteSetOf11HumanAg_Compact_BirthYears.20201120.csv',1,0);


% 16 columns are
%
% 1-5 :  TITER_H1_1918 , TITER_H1_1977 , TITER_H1_1999 , TITER_H1_2007 , TITER_H1_2009 , 
% 6-11:  TITER_H3_1968 , TITER_H3_2003 , TITER_H3_2005 , TITER_H3_2007 , TITER_H3_2009 , TITER_H3_2011 , 
% 12-16: AGE , BIRTHYEAR , LOCATION_CODE , YEAR , 02FL_SAMPLE_ID  

% log-transform the titers
C = log2( B(:,1:11) / 10 );


%%


%
%
[COEFF,SCORE,latent,tsquared,explained,mu] = pca(C, 'Centered', false); % it doesn't seem to matter (much) whether this is true or false

C = log2( B(:,6:11) / 10 );

[COEFF_H3,SCORE_H3,latent_H3,tsquared_H3,explained_H3,mu_H3] = pca(C, 'Centered', false); 

C = log2( B(:,1:5) / 10 );

[COEFF_H1,SCORE_H1,latent_H1,tsquared_H1,explained_H1,mu_H1] = pca(C, 'Centered', false); 

save mfb_workspace_logPC1_analysis.mat

%%

%load mfb_workspace_logPC1_analysis.mat




%%


% %
% %
% %    
% %
% %

aa=log2(B(:,12));
[qq,ii]=sort(aa);
aa_sorted=aa(ii);

plot( aa, SCORE(:,1), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
%set(gca,'XTick', 0:6 );
%set(gca,'XTickLabel', {'1','2','4','8','16','32','64'} );
set(gca,'XTick', [log2(1) log2(2) log2(3) log2(4) log2(5) log2(10) log2(20) log2(40) log2(80) ] );
set(gca,'XTickLabel', {'1','2','3','4','5','10','20','40','80'} );


%yy1 = smooth(aa,SCORE(:,1),0.2,'rloess'); % THIS IS ROBUST LOESS WHICH REMOVES OUTLIERS
%plot(aa,yy1,'.', 'Color', [1.0, 228/255, 225/255]);

sc = SCORE(:,1);
sc_sorted = sc(ii);

yy1 = smooth(aa_sorted,sc_sorted,0.5,'loess');
plot(aa_sorted,yy1,'r-','LineWidth',2);

axis([-0.5 6.8 -14 14])
xlabel('AGE');
ylabel('PC1');

%%

mn = min( SCORE_H3(:,1) );
plot( B(:,12), log( SCORE(:,1) - mn + 0.001), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
xlabel('AGE');
ylabel('LOG PC1');
axis([0 10 -3 5])


%%

% Do the data fitting here

PC1 = SCORE_H3(:,2);
Age = B(:,12);
NMOUTPUT=[];


            options = optimset('Display','off','TolFun',1e-8,'MaxIter', 10000, 'MaxFunEvals',10000);
            iter = 5;

            LB      = [0 0 0 0]; % these 4 parameters in order are H, K, sd, lambda
            UB      = [30 30 30 5];            
            recordPrms = [];
            record_llh = 10000000000;


            for it = 1:iter
                ini = random('unif',LB,UB);
                %LogLikelihoodPC1(ini)
                while isnan(LogLikelihoodPC1(ini))
                   ini = random('unif',LB,UB);
                end

                [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
                NMOUTPUT = [NMOUTPUT; x(1), x(2), x(3), x(4), fval];
                if fval < record_llh
                   record_llh = fval;
                   recordPrms = x;
                end
            end
            
            AR = 1 - exp( - recordPrms(4) );

















%%

% note that below you can get the approximate basis change with
column_means = mean(C);

% make a full matrix of column means
zzz=ones(24402,1)*column_means;

% subtract those column means from the original matrix C
D=C-zzz;

% take some random row in the data set
i=345;

% this should be zero
% and it gives you the transformation from SCORE to your original data
D(i,:)*COEFF - SCORE(i,:);

X = COEFF^(-1);

% to reconstruct the data, do [and this works whether 'Centered' is set to true or false in the PCA]
SCORE(4,:)*X + column_means

% and this will be equal to C(4,:)














