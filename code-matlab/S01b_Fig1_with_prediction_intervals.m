clear all;


% 16 columns are
%
% 1-5 :  TITER_H1_1918 , TITER_H1_1977 , TITER_H1_1999 , TITER_H1_2007 , TITER_H1_2009 , 
% 6-11:  TITER_H3_1968 , TITER_H3_2003 , TITER_H3_2005 , TITER_H3_2007 , TITER_H3_2009 , TITER_H3_2011 , 
% 12-16: AGE , BIRTHYEAR , LOCATION_CODE , YEAR , 02FL_SAMPLE_ID  

% log-transform the titers
%C = log2( B(:,1:11) / 10 );


%%


%
%
%[COEFF,SCORE,latent,tsquared,explained,mu] = pca(C, 'Centered', true); %
% there *are* some minor differences between Centered=true and
% Centered=false; be careful!

%save mfb_workspace_titermatrix.mat

%%

load mfb_workspace_titermatrix.mat


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


lightgreen=[0.92 1.0 0.85];
mediumgreen=[0.3 0.9 0.9];


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
aa=log2(B(:,12)); % this is all the ages
[qq,ii]=sort(aa);
aa_sorted=aa(ii); % sorted ages

plot( aa, SCORE(:,1), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
%set(gca,'XTick', 0:6 );
%set(gca,'XTickLabel', {'1','2','4','8','16','32','64'} );
set(gca,'XTick', [log2(1) log2(2) log2(3) log2(4) log2(5) log2(10) log2(20) log2(40) log2(80) ] );
set(gca,'XTickLabel', {'1','2','3','4','5','10','20','40','80'} );

plot( log2(age_means(:,1)), age_means(:,2), 'bo', 'MarkerFaceColor', 'b' );

%yy1 = smooth(aa,SCORE(:,1),0.2,'rloess'); % THIS IS ROBUST LOESS WHICH REMOVES OUTLIERS
%plot(aa,yy1,'.', 'Color', [1.0, 228/255, 225/255]);

% sort all of the PC1 scores by age, smooth them, and plot them
sc = SCORE(:,1);
sc_sorted = sc(ii);
yy1 = smooth(aa_sorted,sc_sorted,0.5,'loess');
plot(aa_sorted,yy1,'r-','LineWidth',2);


% construct normal prediction intervals by doing local regressions
% in windows of size 0.1 (first arg below) on a log2-age-scale
% last argument gets you the quantile of the prediction interval, so if you
% do 0.95, you will also get the 0.05 quantile, and you are making a 90%
% prediction interval
PRED = get_sigmas_in_windows(  -1.0, 6.8, 0.1 , aa_sorted, sc_sorted, yy1 , 0.90 );

xx=transpose( PRED(:,1) );
yytop = transpose( PRED(:,2)+PRED(:,4) );
yybot = transpose( PRED(:,2)-PRED(:,4) );
h=fill( [xx fliplr(xx)], [yytop fliplr(yybot)], mediumgreen); 
set(h,'facealpha',.16)

% do this to overwrite the unnatura boundaries of the fill
plot( PRED(:,1), PRED(:,2)+PRED(:,4), 'Color', 1-0.5*(1-mediumgreen) );
plot( PRED(:,1), PRED(:,2)-PRED(:,4), 'Color', 1-0.5*(1-mediumgreen) );


axis([-0.5 6.7 -5 22])
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


% construct normal prediction intervals by doing local regressions
% in windows of size 0.1 (first arg below) on a log2-age-scale
% last argument gets you the quantile of the prediction interval, so if you
% do 0.95, you will also get the 0.05 quantile, and you are making a 90%
% prediction interval
PRED = get_sigmas_in_windows(  1916, 2015, 1.0 , bd_sorted, sc_sorted, yy1 , 0.90 );

xx=transpose( PRED(:,1) );
yytop = transpose( PRED(:,2)+PRED(:,4) );
yybot = transpose( PRED(:,2)-PRED(:,4) );
h=fill( [xx fliplr(xx)], [yytop fliplr(yybot)], mediumgreen); 
set(h,'facealpha',.16)

% do this to overwrite the unnatural boundaries of the fill
plot( PRED(:,1), PRED(:,2)+PRED(:,4), 'Color', 1-0.5*(1-mediumgreen) );
plot( PRED(:,1), PRED(:,2)-PRED(:,4), 'Color', 1-0.5*(1-mediumgreen) );

plot( [PRED(1,1) PRED(1,1)], [PRED(1,2)-PRED(1,4),PRED(1,2)+PRED(1,4) ], 'Color', 1-0.5*(1-mediumgreen) );
lst=size(PRED,1);
plot( [PRED(lst,1) PRED(lst,1)], [PRED(lst,2)-PRED(lst,4),PRED(lst,2)+PRED(lst,4) ], 'Color', 1-0.5*(1-mediumgreen) );

plot([1918 1918],[-8 8], 'k' );
plot([1957 1957],[-8 8], 'k' );
plot([1968 1968],[-8 8], 'k' );
plot([1977 1977],[-8 8], 'k' );
plot([2009 2009],[-8 8], 'k' );


axis([1915 2015 -5 5])
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

% construct normal prediction intervals by doing local regressions
% in windows of size 0.1 (first arg below) on a log2-age-scale
% last argument gets you the quantile of the prediction interval, so if you
% do 0.95, you will also get the 0.05 quantile, and you are making a 90%
% prediction interval
PRED = get_sigmas_in_windows(  1916, 2015, 1.0 , bd_sorted, sc_sorted, yy1 , 0.90 );

xx=transpose( PRED(:,1) );
yytop = transpose( PRED(:,2)+PRED(:,4) );
yybot = transpose( PRED(:,2)-PRED(:,4) );
h=fill( [xx fliplr(xx)], [yytop fliplr(yybot)], mediumgreen); 
set(h,'facealpha',.16)

% do this to overwrite the unnatural boundaries of the fill
plot( PRED(:,1), PRED(:,2)+PRED(:,4), 'Color', 1-0.5*(1-mediumgreen) );
plot( PRED(:,1), PRED(:,2)-PRED(:,4), 'Color', 1-0.5*(1-mediumgreen) );

plot( [PRED(1,1) PRED(1,1)], [PRED(1,2)-PRED(1,4),PRED(1,2)+PRED(1,4) ], 'Color', 1-0.5*(1-mediumgreen) );
lst=size(PRED,1);
plot( [PRED(lst,1) PRED(lst,1)], [PRED(lst,2)-PRED(lst,4),PRED(lst,2)+PRED(lst,4) ], 'Color', 1-0.5*(1-mediumgreen) );

plot([1918 1918],[-8 8], 'k' );
plot([1957 1957],[-8 8], 'k' );
plot([1968 1968],[-8 8], 'k' );
plot([1977 1977],[-8 8], 'k' );
plot([2009 2009],[-8 8], 'k' );

axis([1915 2015 -4 4])
xlabel('BIRTH YEAR');
ylabel('PC3');


