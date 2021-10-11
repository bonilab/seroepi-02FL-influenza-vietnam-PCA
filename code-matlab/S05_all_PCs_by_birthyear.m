clear all;


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

Z=[];

    bd=B(:,13);
    [qq,ii]=sort(bd);
    bd_sorted=bd(ii);

for spi=1:11
    

   sc = SCORE(:,spi);
   sc_sorted = sc(ii);

   yy1 = smooth(bd_sorted,sc_sorted,0.5,'loess');

   Z = [Z, yy1];
    
end

save Z.mat Z


%%

% Same smoothing as above, but we are plotting titers by birth year now

ZC = [];

for spi=1:11
    

    cc = C(:,spi);
    cc_sorted = cc(ii);

    yy2 = smooth(bd_sorted,cc_sorted,0.5,'loess');

    ZC = [ZC, yy2];
    
end

save ZC.mat ZC





    

%%

load Z;


lbls={'H1 1918', 'H1 1977', 'H1 1999', 'H1 2007', 'H1 2009', 'H3 1968', 'H3 2003', 'H3 2005', 'H3 2007', 'H3 2009', 'H3 2011'};

ylim=[ -6 2;
       -0.5 0.6;
       -2 1.5;
       -1 1.5;
       -1.5 0.5;
       -0.35 0.3;
       -0.3 0.1;
       -0.15 0.6;
       -0.1 0.2;
       -0.1 0.1;
       -0.05 0.1];


for spi=1:11
    
    subplot(3,4,spi)

    % you don't need bd_sorted below!
    plot( bd, SCORE(:,spi), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
    %set(gca,'XTick', 0:6 );
    %set(gca,'XTickLabel', {'1','2','4','8','16','32','64'} );
    set(gca,'XTick', 1920:15:2020 );
   
    %sc = SCORE(:,spi);
    %sc_sorted = sc(ii);

    %yy1 = smooth(bd_sorted,sc_sorted,0.5,'loess');
    plot(bd_sorted,Z(:,spi),'r-','LineWidth',2);
    
    plot([1918 1918],[-7 2], 'k' );
    plot([1957 1957],[-7 2], 'k' );
    plot([1968 1968],[-7 2], 'k' );
    plot([1977 1977],[-7 2], 'k' );
    plot([2009 2009],[-7 2], 'k' );


    axis([1915 2015 ylim(spi,1) ylim(spi,2)])
    xlabel('BIRTH YEAR');

    ylabel(sprintf('PC %d',spi));

end

% subplot(3,4,12)
% 
% for j=6:6
%    plot(bd_sorted,Z(:,j)); hold on; 
%    axis([1915 2015 -0.2 0.2])
% end







%%



nrz = size(Z,1);
index_of_pcmax = zeros(nrz,1);
index_of_pcmax_excl1 = zeros(nrz,1);

for r =1:nrz
    row1 = Z(r,:);
    [mx,i] = max( row1 );
    index_of_pcmax(r) = i;
    index_of_pcmax_excl1(r) = i;
    
    if i == 1
       row1(i) = -9999;
       [mx,i] = max( row1 );
       index_of_pcmax_excl1(r) = i;
    end
    
end

index_of_pcmin = zeros(nrz,1);
index_of_pcmin_excl = zeros(nrz,1);

for r =1:nrz
    row1 = Z(r,:);
    [mn,i] = min( row1 );
    index_of_pcmin(r) = i;
    index_of_pcmin_excl1(r) = i;

    if i == 1
       row1(i) = 9999;
       [mn,i] = min( row1 );
       index_of_pcmin_excl1(r) = i;
    end
    
end




%%

load ZC;




%%

figure;

% %
% %
% %  SUPP Figure here.  Export as 18" by 10", fixed fonts at 8, eps
% %
% %

for spi=1:11
    
    subplot(3,4,spi)

    plot( bd, C(:,spi), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
    %set(gca,'XTick', 0:6 );
    %set(gca,'XTickLabel', {'1','2','4','8','16','32','64'} );
    set(gca,'XTick', 1920:15:2020 );
   
    %sc = SCORE(:,spi);
    %sc_sorted = sc(ii);

    %yy1 = smooth(bd_sorted,sc_sorted,0.5,'loess');
    plot(bd_sorted,ZC(:,spi),'r-','LineWidth',2);
    
    plot([1918 1918],[-7 8], 'k' );
    plot([1957 1957],[-7 8], 'k' );
    plot([1968 1968],[-7 8], 'k' );
    plot([1977 1977],[-7 8], 'k' );
    plot([2009 2009],[-7 8], 'k' );


    axis([1915 2015 -1 8])
    xlabel('BIRTH YEAR');

    ylabel(sprintf('%s',lbls{spi}));

end




nrz = size(Z,1);
index_of_agmax = zeros(nrz,1);
%index_of_agmax_excl1 = zeros(nrz,1);

for r =1:nrz
    row1 = ZC(r,:);
    [mx,i] = max( row1 );
    index_of_agmax(r) = i;    
end

index_of_agmin = zeros(nrz,1);
%index_of_pcmin_excl = zeros(nrz,1);

for r =1:nrz
    row1 = ZC(r,:);
    [mn,i] = min( row1 );
    index_of_agmin(r) = i;    
end


figure;

subplot(2,2,1)
plot(bd_sorted,index_of_pcmax,'o'); hold on;
plot(bd_sorted,index_of_pcmin,'ro');
axis([1915 2015 0 12])
xlabel('BIRTH YEAR');
ylabel('Index of max PC (blue); Index of min PC (red);');
title('Including PC1');

subplot(2,2,2)
plot(bd_sorted,index_of_pcmax_excl1,'o'); hold on;
plot(bd_sorted,index_of_pcmin_excl1,'ro');
axis([1915 2015 0 12])
xlabel('BIRTH YEAR');
ylabel('Index of max PC (blue); Index of min PC (red);');
title('Excluding PC1');
% subplot(2,2,3)
% plot(bd_sorted,index_of_pcmax,'o'); hold on;
% plot(bd_sorted,index_of_pcmin,'ro');
% axis([1915 2015 0 12])
% xlabel('BIRTH YEAR');
% ylabel('Index of max PC (blue); Index of min PC (red);');
% title('Including PC1');

subplot(2,2,4)
plot(bd_sorted,index_of_agmax,'o'); hold on;
plot(bd_sorted,index_of_agmin,'ro');
axis([1915 2015 0 12])
xlabel('BIRTH YEAR');
set(gca,'YTick', 1:11);
set(gca,'YTickLabel', lbls);
ylabel('Index of max AG (blue); Index of min AG (red);');
title('All Antigens');

%%

% %
% %
% %  FIG 5
% %
% %  Export as 14" by 7",font size 11, eps
% %
% %


figure;
plot(bd_sorted,index_of_pcmax_excl1,'ks','MarkerFaceColor', 'k', 'MarkerSize', 10); hold on;
plot(bd_sorted,index_of_pcmin_excl1,'s','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerSize', 10);
axis([1915 2015 0 12])
xlabel('BIRTH YEAR');
ylabel('Index of max PC (black); Index of min PC (gray);');
%title('Excluding PC1');

% find jump positions in index_of_pcmax1_excl1

J=[];

for r=2:nrz

    if index_of_pcmax_excl1(r) ~= index_of_pcmax_excl1(r-1)
        J = [J, r-1];
    end
end

for r=2:nrz

    if index_of_pcmin_excl1(r) ~= index_of_pcmin_excl1(r-1)
        J = [J, r-1];
    end
end

nb=size(J,2);
for i =1:nb

    plot([bd_sorted(J(i)) bd_sorted(J(i))], [0 12], 'k:' ); 
    
end

text(1930,10,'Exposure to H1-18, H1-77, H1-09, H3-68');
text(1960,10,'Exposure to H1-77, H3-68');

text(1976,8,'H3-03, H3-05');
text(1975,7.5,'H1-77, H1-07, H1-99');

text(1985,10,'H3-03, H3-05');

text(2000,8,'H3');
text(2000,7.5,'H1-09');

text(2004.5,10,'H3-09, H3-11');
text(2005,9.5,'H1-09');

text(2013,2.5,'H3-68');
text(2013,2,'H1-09');

%%

%     bd=B(:,13);
%     [qq,ii]=sort(bd);
%     bd_sorted=bd(ii);
%
% From above, we have bd_sorted already containing all of the sorted birth
% years, and they are sorted by the indicies ii

% now we do a scatter plot of the max antigen without smoothing

% first, sort the matrix B by birth year

B_byr_sorted = B(ii,:);

nrB = size(B_byr_sorted,1);
index_of_maxag_unsmoothed = zeros(nrB,1);
%index_of_pcmax_excl1_unsmoothed = zeros(nrB,1);
index_of_minag_unsmoothed = zeros(nrB,1);
%index_of_pcmin_excl1_unsmoothed = zeros(nrB,1);

for r =1:nrB
    row1 = B_byr_sorted(r,1:11);
    [mx,i] = max( row1 );
    index_of_maxag_unsmoothed(r) = i;
    
end

for r =1:nrB

    row1 = B_byr_sorted(r,1:11);
    [mn,i] = min( row1 );
    index_of_minag_unsmoothed(r) = i;
    
end

%figure;
%plot( bd_sorted , index_of_pcmax_excl1_unsmoothed, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 0.5 );
%axis([1915 2015 0 12]);

%%

% %
% %
% %
% %         RUN THIS BLOCK FOR FIGURE 5 of the main text
% %         export as 15" by 5", fixed fonts at 11, eps
% %
% %


figure;

subplot(1,3,1)
w=0.67; % width is 0.67 years
w1=1915;
w2=1915 + w;
a=0.005;

    plot([1918 1918],[-7 13], 'Color', [1.0 0.5 0.0] ); hold on;
    plot([1957 1957],[-7 13], 'Color', [0.0 0.7 0.0] );
    plot([1968 1968],[-7 13], 'Color', [0.7 0.0 0.0] );
    plot([1977 1977],[-7 13], 'Color', [1.0 0.5 0.0] );
    plot([2009 2009],[-7 13], 'Color', [1.0 0.5 0.0] );

myBlue = [0.1 0.6 1.0];

while w1 < 2015
   
    ind = B_byr_sorted(:,13) >= w1 & B_byr_sorted(:,13) < w2;
    D = B_byr_sorted(ind,:);
    
    x = index_of_maxag_unsmoothed(ind);
    
    size(x);
    
    % should be 1 to 11
    for ag=1:11
    
        y = 0.5*sum(x==ag);
    
        plot([(w1+w2)/2, (w1+w2)/2], [(ag-a*y), (ag+a*y)], 'Color', myBlue, 'LineWidth', 4); hold on;
    
    end
        
    w1 = w1 + w;
    w2 = w2 + w;
    
end



axis([1910 2020 0 12]);
set(gca,'XTick', 1920:20:2000 );
set(gca,'YTick', 1:11);
set(gca,'YTickLabel', lbls);
xlabel('BIRTH YEAR');
ylabel('HIGHEST TITER');


%subplot(1,2,2)
w=0.67; % width is 0.67 years
w1=1915;
w2=1915 + w;
a=0.005;

%     plot([1918 1918],[-7 13], 'Color', [1.0 0.5 0.0] ); hold on;
%     plot([1957 1957],[-7 13], 'Color', [0.0 0.7 0.0] );
%     plot([1968 1968],[-7 13], 'Color', [0.7 0.0 0.0] );
%     plot([1977 1977],[-7 13], 'Color', [1.0 0.5 0.0] );
%     plot([2009 2009],[-7 13], 'Color', [1.0 0.5 0.0] );


while w1 < 2015
   
    ind = B_byr_sorted(:,13) >= w1 & B_byr_sorted(:,13) < w2;
    D = B_byr_sorted(ind,:);
    
    x = index_of_minag_unsmoothed(ind);
    
    size(x);
    
    % should be 1 to 11
    for ag=1:11
    
        y = 0.5*sum(x==ag);
    
        %plot([(w1+w2)/2, (w1+w2)/2], [(ag-a*y), (ag+a*y)], 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
    
    end
        
    w1 = w1 + w;
    w2 = w2 + w;
    
end

%axis([1910 2020 0 12]);
%set(gca,'XTick', 1920:20:2000 );
%xlabel('BIRTH YEAR');
%ylabel('LOWEST TITER');





% first, sort the matrix SCORE by birth year

%     bd=B(:,13);
%     [qq,ii]=sort(bd);
%     bd_sorted=bd(ii);
%
% From above, we have bd_sorted already containing all of the sorted birth
% years, and they are sorted by the indicies ii

% now we do a scatter plot of max PC component without smoothing

% first, sort the matrix B by birth year

SCORE_byr_sorted = SCORE(ii,:);

nrB = size(SCORE_byr_sorted,1);
index_of_pcmax_unsmoothed = zeros(nrB,1);
index_of_pcmax_excl1_unsmoothed = zeros(nrB,1);
index_of_pcmin_unsmoothed = zeros(nrB,1);
index_of_pcmin_excl1_unsmoothed = zeros(nrB,1);

for r =1:nrB
    row1 = SCORE_byr_sorted(r,1:11);
    [mx,i] = max( row1 );
    index_of_pcmax_unsmoothed(r) = i;
    index_of_pcmax_excl1_unsmoothed(r) = i;
    
    if i == 1
       row1(i) = -9999;
       [mx,i] = max( row1 );
       index_of_pcmax_excl1_unsmoothed(r) = i;
    end
    
end

for r =1:nrB
    row1 = SCORE_byr_sorted(r,1:11);
    [mn,i] = min( row1 );
    index_of_pcmin_unsmoothed(r) = i;
    index_of_pcmin_excl1_unsmoothed(r) = i;
    
    if i == 1
       row1(i) = 9999;
       [mn,i] = min( row1 );
       index_of_pcmin_excl1_unsmoothed(r) = i;
    end
    
end

%figure;
%plot( bd_sorted , index_of_pcmax_excl1_unsmoothed, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 0.5 );
%axis([1915 2015 0 12]);



%figure;

subplot(1,3,2)
w=0.67; % width is 0.67 years
w1=1915;
w2=1915 + w;
a=0.005;

    plot([1918 1918],[-7 13], 'Color', [1.0 0.5 0.0] ); hold on;
    plot([1957 1957],[-7 13], 'Color', [0.0 0.7 0.0] );
    plot([1968 1968],[-7 13], 'Color', [0.7 0.0 0.0] );
    plot([1977 1977],[-7 13], 'Color', [1.0 0.5 0.0] );
    plot([2009 2009],[-7 13], 'Color', [1.0 0.5 0.0] );


while w1 < 2015
   
    % STILL NEED TO USE THE B MATRIX BELOW BC SCORE DOES NOT HAVE BIRTH
    % YEARS STORED IN IT
    ind = B_byr_sorted(:,13) >= w1 & B_byr_sorted(:,13) < w2;
    %D = B_byr_sorted(ind,:);
    
    x = index_of_pcmax_excl1_unsmoothed(ind);
    
    size(x);
    
    % should be 2 to 11
    for pc=2:11
    
        y = 0.5*sum(x==pc);
    
        plot([(w1+w2)/2, (w1+w2)/2], [(pc-a*y), (pc+a*y)], 'k', 'LineWidth', 4); hold on;
    
    end
        
    w1 = w1 + w;
    w2 = w2 + w;
    
end



axis([1910 2020 0 12]);
set(gca,'XTick', 1920:20:2000 );
xlabel('BIRTH YEAR');
ylabel({'MAXIMUM POSITIVE','PRINCIPAL COMPONENT'});


subplot(1,3,3)
w=0.67; % width is 0.67 years
w1=1915;
w2=1915 + w;
a=0.005;

    plot([1918 1918],[-7 13], 'Color', [1.0 0.5 0.0] ); hold on;
    plot([1957 1957],[-7 13], 'Color', [0.0 0.7 0.0] );
    plot([1968 1968],[-7 13], 'Color', [0.7 0.0 0.0] );
    plot([1977 1977],[-7 13], 'Color', [1.0 0.5 0.0] );
    plot([2009 2009],[-7 13], 'Color', [1.0 0.5 0.0] );


while w1 < 2015
   
    % STILL NEED TO USE THE B MATRIX BELOW BC SCORE DOES NOT HAVE BIRTH
    % YEARS STORED IN IT
    ind = B_byr_sorted(:,13) >= w1 & B_byr_sorted(:,13) < w2;
    % D = B_byr_sorted(ind,:);
    
    x = index_of_pcmin_excl1_unsmoothed(ind);
    
    size(x);
    
    % should be 2 to 11
    for pc=2:11
    
        y = 0.5*sum(x==pc);
    
        plot([(w1+w2)/2, (w1+w2)/2], [(pc-a*y), (pc+a*y)], 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
    
    end
        
    w1 = w1 + w;
    w2 = w2 + w;
    
end

axis([1910 2020 0 12]);
set(gca,'XTick', 1920:20:2000 );
xlabel('BIRTH YEAR');
ylabel({'MAXIMUM NEGATIVE','PRINCIPAL COMPONENT'});


