clear all;

strBase = 'AR15.Sigma0p88.nb0p3.rep03';    
numfiles=1;        

for i=1:numfiles

        strInfile = ['./B_Simulated/BMatrix.' strBase '.csv'];
        %strOutfile = ['./AgePC1_Simulated/AgePC1.' strBase '.csv'];

        B = csvread( strInfile ,1,0);
        % 16 columns are (except first five are -99, and 13 to 16 are -99
        %
        % 1-5 :  TITER_H1_1918 , TITER_H1_1977 , TITER_H1_1999 , TITER_H1_2007 , TITER_H1_2009 , 
        % 6-11:  TITER_H3_1968 , TITER_H3_2003 , TITER_H3_2005 , TITER_H3_2007 , TITER_H3_2009 , TITER_H3_2011 , 
        % 12-16: AGE , BIRTHYEAR , LOCATION_CODE , YEAR , 02FL_SAMPLE_ID  

        % log-transform the H3N2 titers
        C = log2( B(:,6:11) / 10 );

        [COEFF,SCORE,latent,tsquared,explained,mu] = pca(C, 'Centered', false); % it doesn't seem to matter whether this is true or false

        % compute the PC-averages by age and birthyear
        D = [SCORE, B(:,12)];
        numrows=size(D,1);
        
        % make an empty vector for the PC1 to PC4 means by age
        age_means=ones(90,5);
        age_means = age_means * (-99);
        age_means(:,1) = transpose(1:90);
        
        for aa=1:90
        
            E = D( D(:,7)>=aa-0.5 & D(:,7)<aa+0.5, : );
            col_means = mean(E);
            age_means(aa,2:5) = col_means(1:4);
            
        end

        aa=log2(B(:,12));
        [qq,ii]=sort(aa);
        aa_sorted=aa(ii);
        
        minPC1 = min( SCORE(:,1) );

        plot( aa, SCORE(:,1), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
        plot( [-2 8], [minPC1 minPC1], 'k-' );
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

        %axis([-0.5 6.8 (minPC1-0.2) 14])
        %axis([-0.5 6.8 7 11])
        xlabel('AGE');
        ylabel('PC1');

        WM = [B(:,12), SCORE(:,1)];
        %figure;
        %scatter( log2(WM(:,1)), WM(:,2) );
        %writematrix( WM, strOutfile );
end
