clear all;
%A = csvread('./CrossValidation/FINAL_InferredAR.AL.csv',0,0);
A = csvread('./CrossValidation/FINAL_InferredAR.AL.withCI.csv',0,0);


nr=size(A,1);
col_zeros = zeros(nr,1);
A = [A, col_zeros];

% columns
% 1 subset being analyzed
% 2 PC basis from this subset is used
% 3 location (0=combined for all locs)
% 4 attack rate
% 5-6 : lower, upper confidence intervals on attack rate
% 7 subset being analyzed -- use this for plotting as it has some useful
%                            horizontal offsets to improve visibility

my_greenblue_1 = [161,218,180]/256;
my_greenblue_2 = [65,182,196]/256;
my_greenblue_3 = [44,127,184]/256;
my_greenblue_4 = [37,52,148]/256;
my_greenblue_5 = [20,22,118]/256; % made up by MFB


% column 5 is the subset number .... but it is offset left/right just a bit
% for better visibility
for r = 1:nr
   
    x = (A(r,2)-3) * 0.1;
    A(r,7) = A(r,1) + x; 
    
end

% draw the background horizontal lines showing the correct AR
for k = 1:5
    
    D = A( A(:,1)==k & A(:,2)==k , : );
    %size(D)
    plot([k-0.4,k+0.4],[D(4) D(4)], '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on; 
    
end

% draw the vertical confidence interval lines -- su=subset, ba=basis
for su=1:5
    for ba=1:5
        
        D = A( A(:,1)==su & A(:,2)==ba , : );
        %size(D)
        plot([D(7) D(7)],[D(5) D(6)], 'k-' ); hold on;
        
        
    end
end


scatter( A(:,7), A(:,4), 'MarkerFaceColor', my_greenblue_2, 'MarkerEdgeColor', my_greenblue_5, 'SizeData', 120, 'MarkerFaceAlpha', 0.5 ); hold on;

axis([-0.6 6 0.20 .35]);
set(gca, 'XTick', 1:5 );
set(gca, 'XTickLabel', {'SUBSET 1', 'SUBSET 2', 'SUBSET 3', 'SUBSET 4', 'SUBSET 5'} );
grid on;
ylabel({'INFERRED ATTACK RATE IN','CROSS-VALIDATION EXERCISE'});



