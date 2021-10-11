clear all;

B = csvread('../data/TITERS_02FL_CompleteSetOf11HumanAg_Compact_BirthYears.20201120.csv',1,0);


% 16 columns in the above file are
%
% 1-5 :  TITER_H1_1918 , TITER_H1_1977 , TITER_H1_1999 , TITER_H1_2007 , TITER_H1_2009 , 
% 6-11:  TITER_H3_1968 , TITER_H3_2003 , TITER_H3_2005 , TITER_H3_2007 , TITER_H3_2009 , TITER_H3_2011 , 
% 12-16: AGE , BIRTHYEAR , LOCATION_CODE , YEAR , 02FL_SAMPLE_ID  

% log-transform the titers
C = log2( B(:,1:11) / 10 );


%%


%
%
[COEFF,SCORE,latent,tsquared,explained,mu] = pca(C, 'Centered', false); % it doesn't seem to matter whether this is true or false

save mfb_workspace_titermatrix.mat

%%

%load mfb_workspace_titermatrix.mat


%%

% compute the PC-averages by age and birthyear

D = [SCORE, B(: ,12:13)];
numrows=size(D,1);

% just do the first four PCs
age_means=ones(90,5);
age_means = age_means * (-99);
age_means(:,1) = transpose(1:90);

for aa=1:90

    E = D( D(:,12)>=aa-0.5 & D(:,12)<aa+0.5, : );
    col_means = mean(E);
    age_means(aa,2:5) = col_means(1:4);
    
end

% just do the first four PCs
byear_means=ones(99,5);
byear_means = byear_means * (-99);
byear_means(:,1) = transpose(1916:2014);

for by=1916:2014

    E = D( D(:,13)>=by & D(:,13)<by+1, : );
    col_means = mean(E);
    index=by-1915;
    byear_means(index,2:5) = col_means(1:4);
    
end






%%


% %
% %
% %    THIS CODE BLOCK MAKES FIGURE 2
% %
% %

% %     FIRST, MANUALLY MOVE ANTIGEN LABELS SO THEY ARE VISIBLE
% %     EXPORT AS 20" BY 12", FIXED FONTS AT 10.  APPLY TO FIGURE.  EXPORT
% %     AS EPS

lbls={'H1 1918', 'H1 1977', 'H1 1999', 'H1 2007', 'H1 2009', 'H3 1968', 'H3 2003', 'H3 2005', 'H3 2007', 'H3 2009', 'H3 2011'};

subplot(2,3,1)
biplot( COEFF(:,1:2), 'Scores', SCORE(1,1:2), 'VarLabels', lbls );
xlabel( sprintf('PC1 (%1.1f%% of total variance)', 100*latent(1)/sum(latent) ) );
ylabel( sprintf('PC2 (%1.1f%% of total variance)', 100*latent(2)/sum(latent) ) );
subplot(2,3,2)
biplot( COEFF(:,2:3), 'Scores', SCORE(1,2:3), 'VarLabels', lbls );
xlabel( sprintf('PC2 (%1.1f%% of total variance)', 100*latent(2)/sum(latent) ) );
ylabel( sprintf('PC3 (%1.1f%% of total variance)', 100*latent(3)/sum(latent) ) );
subplot(2,3,3)
biplot( COEFF(:,3:4), 'Scores', SCORE(1,3:4), 'VarLabels', lbls );
xlabel( sprintf('PC3 (%1.1f%% of total variance)', 100*latent(3)/sum(latent) ) );
ylabel( sprintf('PC4 (%1.1f%% of total variance)', 100*latent(4)/sum(latent) ) );



subplot(2,3,4)
aa=log2(B(:,12));
[qq,ii]=sort(aa);
aa_sorted=aa(ii);

plot( aa, SCORE(:,1), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
%set(gca,'XTick', 0:6 );
%set(gca,'XTickLabel', {'1','2','4','8','16','32','64'} );
set(gca,'XTick', [log2(1) log2(2) log2(3) log2(4) log2(5) log2(10) log2(20) log2(40) log2(80) ] );
set(gca,'XTickLabel', {'1','2','3','4','5','10','20','40','80'} );

plot( log2(age_means(:,1)), age_means(:,2), 'bo', 'MarkerFaceColor', 'b' );

%yy1 = smooth(aa,SCORE(:,1),0.2,'rloess'); % THIS IS ROBUST LOESS WHICH REMOVES OUTLIERS
%plot(aa,yy1,'.', 'Color', [1.0, 228/255, 225/255]);

sc = SCORE(:,1);
sc_sorted = sc(ii);

yy1 = smooth(aa_sorted,sc_sorted,0.5,'loess');
plot(aa_sorted,yy1,'r-','LineWidth',2);

axis([-0.5 6.8 -14 14])
xlabel('AGE');
ylabel('PC1');



subplot(2,3,5)
bd=B(:,13);
[qq,ii]=sort(bd);
bd_sorted=bd(ii);

plot( bd, SCORE(:,2), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
set(gca,'XTick', 1920:15:2020 );
%set(gca,'XTickLabel', {'1','2','3','4','5','10','20','30','50','80'} );

plot( byear_means(:,1), byear_means(:,3), 'bo', 'MarkerFaceColor', 'b' );

%yy1 = smooth(bd,SCORE(:,2),0.2,'rloess'); % THIS IS ROBUST LOESS WHICH REMOVES OUTLIERS
%plot(bd,yy1,'.', 'Color', [1.0, 228/255, 225/255]);

sc = SCORE(:,2);
sc_sorted = sc(ii);

yy1 = smooth(bd_sorted,sc_sorted,0.5,'loess');
plot(bd_sorted,yy1,'r-','LineWidth',2);

plot([1918 1918],[-2 2], 'k' );
plot([1957 1957],[-2 2], 'k' );
plot([1968 1968],[-2 2], 'k' );
plot([1977 1977],[-2 2], 'k' );
plot([2009 2009],[-2 2], 'k' );


axis([1915 2015 -2 2])
xlabel('BIRTH YEAR');
ylabel('PC2');

subplot(2,3,6)

plot( bd, SCORE(:,3), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
set(gca,'XTick', 1920:15:2020 );

plot( byear_means(:,1), byear_means(:,4), 'bo', 'MarkerFaceColor', 'b' );


%yy1 = smooth(bd,SCORE(:,3),0.2,'rloess'); % THIS IS ROBUST LOESS WHICH REMOVES OUTLIERS
%plot(bd,yy1,'.', 'Color', [1.0, 228/255, 225/255]);

sc = SCORE(:,3);
sc_sorted = sc(ii);

yy1 = smooth(bd_sorted,sc_sorted,0.5,'loess');
plot(bd_sorted,yy1,'r-','LineWidth',2);

plot([1918 1918],[-2 2], 'k' );
plot([1957 1957],[-2 2], 'k' );
plot([1968 1968],[-2 2], 'k' );
plot([1977 1977],[-2 2], 'k' );
plot([2009 2009],[-2 2], 'k' );

axis([1915 2015 -2 2])
xlabel('BIRTH YEAR');
ylabel('PC3');


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

%%
figure;

for j=1:11

    a=-column_means(j);
    b=COEFF(j,1)*(7.5-column_means(j));
    plot([0 7.5],[a b], 'k-'); hold on;
    
    
end


%%

figure;

for spi=1:10

    subplot(2,5,spi)
   
    biplot( COEFF(:,spi:(spi+1)), 'Scores', SCORE(1,spi:(spi+1)), 'VarLabels', lbls );
    
    xlabel(sprintf('%d',spi));
    ylabel(sprintf('%d',spi+1));
    
end

%%
figure;
j=3;
subplot(1,2,1)
biplot( COEFF(:,(j-1):j), 'Scores', SCORE(1,(j-1):j), 'VarLabels', lbls );
    xlabel(sprintf('%d',j-1));
    ylabel(sprintf('%d',j));
subplot(1,2,2)
biplot( COEFF(:,j:(j+1)), 'Scores', SCORE(1,j:(j+1)), 'VarLabels', lbls );
    xlabel(sprintf('%d',j));
    ylabel(sprintf('%d',j+1));



%%

figure;
plot( B(:,12), B(:,6), 'o' );

BB06 = B(:,6);
BB12 = B(:,12);

A50 = B( B(:,12) > 50.0 & B(:,12) < 55.0 , : );
A10 = B( B(:,12) > 10.0 & B(:,12) < 15.0 , : );

A40 = B( B(:,12) > 40.0 & B(:,12) < 45.0 , : );
A30 = B( B(:,12) > 30.0 & B(:,12) < 35.0 , : );
A20 = B( B(:,12) > 20.0 & B(:,12) < 25.0 , : );
mean( A40(:,6) )





















