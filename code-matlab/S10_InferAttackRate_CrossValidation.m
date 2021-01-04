%clc;
clear all;
%==============================================
global PC1 Age

all_files = {};


for basis = 1:5
   for subset=1:5

            new_str = sprintf('AgePC1.AL.Subset%d.UsingBasis%d.csv', subset, basis);
            all_files{end+1} = new_str;    
       
   end
end



numfiles=size(all_files,2);  


%strBase = 'AR15.Sigma0p85.nb0p9.rep07';

for i=1:numfiles

            i
            strBase = all_files{i}
            strInfile  = ['./CrossValidation/' strBase ];
            strOutfile = ['./CrossValidation/InferredAR.' strBase ];

            inputSubset = str2double( strBase(17) );
            inputBasis = str2double( strBase(29) );

            data = csvread( strInfile , 0 , 0 );
            PC1 = data(:,2);
            Age = data(:,1);
            %Age = Age - 0.5;
            %===========Search Routine Settings =================
            %options = optimset('Display','iter','TolFun',1e-8,'MaxIter', 10000, 'MaxFunEvals',10000);
            options = optimset('Display','off','TolFun',1e-8,'MaxIter', 10000, 'MaxFunEvals',10000);
            iter = 10;

            LB      = [0 0 0 0];
            UB      = [30 30 30 5];
            % initial = random('unif',LB,UB)

            %====================================================

             recordPrms = [];
             record_llh = 10000000000;


             for it = 1:iter
                 ini = random('unif',LB,UB);
                 while isnan(LogLikelihoodPC1(ini))
                    ini = random('unif',LB,UB);
                 end

                 [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
                 if fval < record_llh
                    record_llh = fval;
                    recordPrms = x;
                 end
             end

            AR = 1 - exp( - recordPrms(4) );

            dlmwrite(strOutfile,[inputSubset inputBasis 0 AR],'delimiter', ',','precision', 6); % that zero means that it is for all sites
 
end



