clear all;

B = csvread('../data/TITERS_02FL_CompleteSetOf11HumanAg_Compact_BirthYears.20201120.csv',1,0);

% 17 columns are
%
% 1-5 :  TITER_H1_1918 , TITER_H1_1977 , TITER_H1_1999 , TITER_H1_2007 , TITER_H1_2009 , 
% 6-11:  TITER_H3_1968 , TITER_H3_2003 , TITER_H3_2005 , TITER_H3_2007 , TITER_H3_2009 , TITER_H3_2011 , 
% 12-16: AGE , BIRTHYEAR , LOCATION_CODE , YEAR , 02FL_SAMPLE_ID  
% 17   : assigned random subgroup, from one to five (this col is added below)

nr = size(B,1);
m = floor(nr/5); % number of samples in each subset
m2=nr-4*m;

% assign each individual into one of five random subgroups
group_vector = [ones(m,1); 2*ones(m,1); 3*ones(m,1); 4*ones(m,1); 5*ones(m2,1)];

group_assignments = group_vector(randperm(nr));
B = [ B , group_assignments ];

A{1} = B( B(:,17)==1 , : );
A{2} = B( B(:,17)==2 , : );
A{3} = B( B(:,17)==3 , : );
A{4} = B( B(:,17)==4 , : );
A{5} = B( B(:,17)==5 , : );

both_subtypes = 1:11;
h1_only = 1:5;
h3_only = 6:11;

selected_columns = h3_only;

C{1} = log2( A{1}( : , selected_columns ) / 10 );
C{2} = log2( A{2}( : , selected_columns ) / 10 );
C{3} = log2( A{3}( : , selected_columns ) / 10 );
C{4} = log2( A{4}( : , selected_columns ) / 10 );
C{5} = log2( A{5}( : , selected_columns ) / 10 );



%%

COEFF = {};
SCORE = {};



for i =1:5
     
    [COEFF{i},SCORE{i},latent,tsquared,explained,mu] = pca(C{i}, 'Centered', true); % unlike previous versions of Matlab, it does matter whether this is 
                                                                                    % set to true or false; make sure it is set to true
end

save mfb_workspace_titermatrix_subsets_H3only.mat

%%

%load mfb_workspace_titermatrix_subsets_H3only.mat




%%

% note that below you can get the approximate basis change with
column_means = mean(C{1});

% make a full matrix of column means ....nr=24402
% below, it's a column vector times a row vector, so the result is a large
% 4880 x 11 matrix
zzz=ones(m,1)*column_means;

% subtract those column means from the original matrix C
% this now holds the centered log-titer values
D=C{1}-zzz;

% take some random row in the data set
i=345;

% this should be zero
% and it gives you the transformation from SCORE to your original data
D(i,:)*COEFF{1} - SCORE{1}(i,:)

X = COEFF{1}^(-1);

% to reconstruct a data point, do this
SCORE{1}(220,:)*X + column_means
C{1}(220,:)
% and this will be equal to C(220,:), i.e. the log-titer values for
% individual #220


%%

% this is the set of basis vectors that we're using;
% in other words, if i=3, that means we do a PCA on subset 3, get the basis
% vectors, and then use those to transform the data in all the other
% subsets (including subset 3)
%
for i = 1:5

    
        % this is the subset we are examining and for which
        % we are computing an attack rate
        for j=1:5

                age_col = A{j}(:,12);
                ind_hc = A{j}(:,14) == 1;
                ind_kh = A{j}(:,14) == 2;
                ind_hu = A{j}(:,14) == 3;

                % number of subset rows
                nsr = size(A{j},1);

                % take the real data from group j and center it
                column_means = mean(C{j});
                D = C{j} - ones(nsr,1)*column_means;

                SCORE_CRV = D * COEFF{i};   % use the basis transformations from coefficient matrix i, because 
                                            % you are generating PC values for
                                            % group j based on the coefficients
                                            % from PCA on group i
                pc1_col = SCORE_CRV(:,1);

                strOutfile = sprintf('./CrossValidation/AgePC1.AL.Subset%d.UsingBasis%d.csv', j, i);                                    
                writematrix( [age_col , pc1_col ], strOutfile );

                strOutfile = sprintf('./CrossValidation/AgePC1.HC.Subset%d.UsingBasis%d.csv', j, i);                                    
                writematrix( [age_col(ind_hc) , pc1_col(ind_hc) ], strOutfile );

                strOutfile = sprintf('./CrossValidation/AgePC1.KH.Subset%d.UsingBasis%d.csv', j, i);                                    
                writematrix( [age_col(ind_kh) , pc1_col(ind_kh) ], strOutfile );

                strOutfile = sprintf('./CrossValidation/AgePC1.HU.Subset%d.UsingBasis%d.csv', j, i);                                    
                writematrix( [age_col(ind_hu) , pc1_col(ind_hu) ], strOutfile );


        end

end






















