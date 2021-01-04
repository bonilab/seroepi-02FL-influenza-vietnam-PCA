%clc;
clear all;
%==============================================
global PC1 Age

all_files = {};

ars=[5 10 15 20 25 30];
sigmas=[0.88];
reps=1:30;
betas = [0.1 0.3 0.5];

for ar=ars
    for sigma=sigmas
        for rep=reps
            for beta=betas
           
                    new_str = sprintf('AR%02d.Sigma0p%02d.nb0p%1d.rep%02d', ar, 100*sigma, beta*10, rep);
                    %new_str = sprintf('AR%02d.Sigma0p%02d.nb1p0.rep%02d', ar, 100*sigma, rep);
                    all_files{end+1} = new_str;    
                    
            end
        end
    end
end

numfiles=size(all_files,2);  


%strBase = 'AR15.Sigma0p85.nb0p9.rep07';

for i=1:numfiles

            i
            strBase = all_files{i}
            strInfile  = ['./AgePC1_Simulated/AgePC1.' strBase '.csv'];
            strOutfile = ['./AgePC1_Simulated/Inferred.' strBase '.csv'];

            inputAR = str2double( strBase(3:4) );
            inputSigma = str2double( strBase(13:14) ) / 100;
            inputBeta = str2double(strBase(18)) + ( str2double(strBase(20))/10 );
            inputBeta = (-1)*inputBeta;
            inputRep = str2double( strBase(25:26) );

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

            dlmwrite(strOutfile,[inputAR inputSigma inputBeta inputRep recordPrms record_llh AR],'delimiter', ',','precision', 6);
 
end



