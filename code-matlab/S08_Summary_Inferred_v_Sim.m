clear all;



% THESE HAVE 10 REPLICATES
%
%A = csvread( './AgePC1_Simulated/FINAL_Inferred.summary.500K.beta135.take3.afterbugfix.csv',0,0);
%A = csvread( './AgePC1_Simulated/FINAL_Inferred.summary.500K.beta135.boost7p26.csv',0,0);



% THESE HAVE 50 REPLICATES
%
A = csvread( './AgePC1_Simulated/FINAL_Inferred.summary.500K.beta135.rep50.csv',0,0);
%A = csvread( './AgePC1_Simulated/FINAL_Inferred.summary.500K.beta135.boost7p26.rep50.take2.csv',0,0);
%




% columns are
%
% 1 the attack rate in the simulation (need to divide by 100 here)
% 2 the sigma value in the simulation
% 3 the beta value in the simulation
% 4 the simulation replicate number
% 5-9 can ignore
% 10 the inferred attack rate

my_greenblue_1 = [161,218,180]/256;
my_greenblue_2 = [65,182,196]/256;
my_greenblue_3 = [44,127,184]/256;
my_greenblue_4 = [37,52,148]/256;
my_greenblue_5 = [20,22,118]/256; % made up by MFB

mycolors={my_greenblue_1, my_greenblue_2, my_greenblue_3, my_greenblue_4, my_greenblue_5};

%sigmas = [0.50 0.60 0.70];
% sigmas = [0.4 0.6 0.8];



subplot(1,2,1)

plot([0.05 0.05],[0 1], 'k:'); hold on;
plot([0.10 0.10],[0 1], 'k:'); 
plot([0.15 0.15],[0 1], 'k:'); 
plot([0.20 0.20],[0 1], 'k:'); 
plot([0.25 0.25],[0 1], 'k:'); 
plot([0.30 0.30],[0 1], 'k:'); 
plot([0.35 0.35],[0 1], 'k:'); 

%axis([0 0.40 0 0.40]);


betas = [-0.5 -0.3 -0.1];
%betas = [-0.3 -0.1];

nb = size(betas,2);
sd_circle_size = 50;

for i=1:nb
    
    %B = A( A(:,2)==sigmas(i), : );
    B = A( A(:,3)==betas(i), : );
    %B = B( B(:,3)==betas(2), : );
    
    %x = (sigmas(i)-0.6)*0.04;
    x = (betas(i)+0.3)*0.05;
    
    plot( [0.0 0.50], [0.0 0.50], 'k-' ); hold on;
    %h=plot( A(:,1)/100, A(:,10), 'ko', 'MarkerFaceColor', my_greenblue_1, 'MarkerEdgeColor', my_greenblue_1, 'MarkerSize', 9 );
    h=scatter( B(:,1)/100 + x, B(:,10), 'MarkerFaceColor', mycolors{i}, 'MarkerEdgeColor', mycolors{i}, 'SizeData', sd_circle_size, 'MarkerFaceAlpha', 0.5 );
    %h=scatter( B(:,1)/100 + x, B(:,10), 'MarkerFaceColor', mycolors{i}, 'MarkerEdgeColor', mycolors{i}, 'SizeData', 40, 'MarkerFaceAlpha', 0.5 );
    %alpha(h,0.5);

    xlabel('SIMULATED ATTACK RATE');
    ylabel({'INFERRED ATTACK RATE','FROM PRINCIPAL COMPONENTS ANALYSIS'});

end
axis([0 0.4 0 0.4]);
%title('beta = -0.1; sigma = 0.4, 0.6, 0.8 from left to right');

subplot(1,2,2)
plot([0.05 0.05],[0 1], 'k:'); hold on;
plot([0.10 0.10],[0 1], 'k:'); 
plot([0.15 0.15],[0 1], 'k:'); 
plot([0.20 0.20],[0 1], 'k:'); 
plot([0.25 0.25],[0 1], 'k:'); 
plot([0.30 0.30],[0 1], 'k:'); 
plot([0.35 0.35],[0 1], 'k:'); 

for i=1:nb
    for ar=5:5:30
    
                B = A( A(:,3)==betas(i), : );
                C = B( B(:,1)==ar, : );
                qs = quantile(C(:,10), [.25 .5 .75]);

                x = (betas(i)+0.3)*0.03;

                plot( [0.0 0.50], [0.0 0.50], 'k-' ); hold on;

                %h=scatter( B(:,1)/100 + x, B(:,10), 'MarkerFaceColor', mycolors{i}, 'MarkerEdgeColor', mycolors{i}, 'SizeData', 110, 'MarkerFaceAlpha', 0.5 );
                plot([(ar/100)+x,(ar/100)+x],[qs(1), qs(3)], 'Color', mycolors{i}, 'LineWidth', 3 ); hold on;
                plot( ar/100 + x, qs(2), 'o', 'MarkerFaceColor', mycolors{i}, 'MarkerEdgeColor', mycolors{i}, 'MarkerSize', 7 );


                %alpha(h,0.5);

                xlabel('SIMULATED ATTACK RATE');
                ylabel({'INFERRED ATTACK RATE','FROM PRINCIPAL COMPONENTS ANALYSIS'});

    end
end

axis([0 0.4 0 0.4]);

% export as 16" by 6" w 12pt font, and to jpg
